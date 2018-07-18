#pragma once

#include "types.hpp"
#include "geometry.hpp"
#include "photometric_consistency.hpp"

static FloatT TriangulationAngle(const SrcView &src, const Vector3 &X)
{
    const Vector3 dX = X - src.C;

    FloatT c = dX.dot(X) / (dX.norm() * X.norm());
    if (c > +1) c = +1;
    if (c < -1) c = -1;

    const FloatT triangulation_angle = std::abs(std::acos(c));

    return triangulation_angle;
}

static FloatT TriangulationPrior(const FloatT triangulation_angle, const FloatT min_triangulation_prior)
{
    const FloatT trinagulation_angle_diff = (std::min(triangulation_angle, min_triangulation_prior) - min_triangulation_prior);

    return 1 - (trinagulation_angle_diff * trinagulation_angle_diff) / (min_triangulation_prior * min_triangulation_prior);
}

static FloatT TriangulationPrior(const SrcView &src, const Vector3 &X, const FloatT min_triangulation_prior)
{
    return TriangulationPrior(TriangulationAngle(src, X), min_triangulation_prior);
}

static FloatT ResolutionPrior(const RefView &ref, const SrcView &src, const Normal &n, const Vector3 &X, const Vector2i &x, int window_size)
{
    const FloatT d = X.dot(n);

    const Normal ray_A = ref.ray(convertToFloat((x + Vector2i(-window_size, -window_size)).eval()));  // top-left
    const Normal ray_B = ref.ray(convertToFloat((x + Vector2i(+window_size, -window_size)).eval()));  // top-right
    const Normal ray_C = ref.ray(convertToFloat((x + Vector2i(+window_size, +window_size)).eval()));  // bottom-right
    const Normal ray_D = ref.ray(convertToFloat((x + Vector2i(-window_size, +window_size)).eval()));  // bottom-left

    const Vector2 x_src_A = src.project((d / ray_A.dot(n)) * ray_A);
    const Vector2 x_src_B = src.project((d / ray_B.dot(n)) * ray_B);
    const Vector2 x_src_C = src.project((d / ray_C.dot(n)) * ray_C);
    const Vector2 x_src_D = src.project((d / ray_D.dot(n)) * ray_D);

    const FloatT ref_area = (2 * window_size + 1) * (2 * window_size + 1);
    const FloatT src_area = QuadrilateralArea(x_src_A, x_src_B, x_src_C, x_src_D);

    if (ref_area < src_area)
        return ref_area / src_area;
    else
        return src_area / ref_area;
}

static FloatT ResolutionPrior(const RefView &ref, const SrcView &src, const Normal &n, const FloatT depth, const Vector2i &x, int window_size)
{
    return ResolutionPrior(ref, src, n, ref.unproject(convertToFloat(x), depth), x, window_size);
}

static const FloatT IncidentAngle(const SrcView &src, const Vector3 &X, const Vector3 &n)
{
    const Vector3 dX = src.C - X;

    FloatT c = dX.dot(n) / (dX.norm() * n.norm());
    if (c > +1) c = +1;
    if (c < -1) c = -1;

    const FloatT incident_angle = std::abs(std::acos(c));

    return incident_angle;
}

static FloatT IncidentPrior(const FloatT incident_angle, FloatT incident_prior_sigma)
{
    return std::exp(- (incident_angle * incident_angle) / (2 * incident_prior_sigma * incident_prior_sigma));
}

static FloatT IncidentPrior(const SrcView &src, const Vector3 &X, const Vector3 &n, FloatT incident_prior_sigma)
{
    return IncidentPrior(IncidentAngle(src, X, n), incident_prior_sigma);
}

// TODO: Make this a class with above functions as members
static FloatT PriorProbability(const RefView &ref, const SrcView &src, const Normal &n, const FloatT depth, const Vector2i &x, int window_size, const FloatT incident_prior_sigma, const FloatT min_triangulation_prior)
{
    const Vector3 X = ref.unproject(convertToFloat(x), depth);

    FloatT triangulation_prior = TriangulationPrior(src, X, min_triangulation_prior);
    FloatT resolution_prior = ResolutionPrior(ref, src, n, X, x, window_size);
    FloatT incident_prior = IncidentPrior(src, X, n, incident_prior_sigma);

    if (std::isnan(triangulation_prior))
    {
        std::cerr << "triangulation_prior: " << triangulation_prior << std::endl;
        triangulation_prior = 0;
    }

    if (std::isnan(resolution_prior))
    {
        std::cerr << "resolution_prior: " << resolution_prior << std::endl;
        resolution_prior = 0;
    }

    if (std::isnan(incident_prior))
    {
        std::cerr << "incident_prior: " << incident_prior << std::endl;
        incident_prior = 0;
    }

    return triangulation_prior * resolution_prior * incident_prior;
}

static FloatT StateTransitionProbability(int current_transition, int total_transitions)
{
//    const FloatT initial_probability = 0.0; /// TODO: colmap's default value

    //const FloatT initial_probability = 0.5; // TOOD: proper paper based value
    const FloatT initial_probability = 0.0; // TOOD: testing value
    return 1 - ((1 - initial_probability) * current_transition / total_transitions + initial_probability);
}

static FloatT SelectionProbability(const FloatT forward_message, const FloatT backward_message, const FloatT previous_probability, const FloatT state_transition_probability)
{
    const FloatT sel_prob_0 = (1 - forward_message) * (1 - backward_message);
    const FloatT sel_prob_1 = forward_message * backward_message;

    const FloatT current_probability = sel_prob_1 / (sel_prob_0 + sel_prob_1); // Eq: 11
    const FloatT selection_probability = (1 - state_transition_probability) * previous_probability + state_transition_probability * current_probability;

    return selection_probability;
}

struct MessageCalculator
{
    MessageCalculator(const FloatT color_similarity_sigma, const FloatT no_change_prob = 0.999)
        : inv_sigma_sq(1 / (2 * color_similarity_sigma * color_similarity_sigma))
        , visibility_probability_normalization(1 / VisibilityProbabilityNormalization(color_similarity_sigma))
        , no_change_prob(no_change_prob)
        , change_prob(1 - no_change_prob)
    {

    }

    FloatT forwardMessage(const FloatT color_similarity_ncc, const FloatT prev_forward_message)
    {
        // const FloatT no_change_prob = 0.99999; /// TODO: colmap's default value

        const FloatT visibility_probability_0 = 0.5f;
        const FloatT visibility_probability_1 = visibilityProbability(color_similarity_ncc);

        const FloatT message_prob_0 = (prev_forward_message * change_prob + (1 - prev_forward_message) * no_change_prob) * visibility_probability_0;
        const FloatT message_prob_1 = (prev_forward_message * no_change_prob + (1 - prev_forward_message) * change_prob) * visibility_probability_1;

        return message_prob_1 / (message_prob_0 + message_prob_1);
    }


    FloatT backwardMessage(const FloatT color_similarity_ncc, const FloatT prev_backward_message)
    {
        // const FloatT no_change_prob = 0.99999; /// TODO: colmap's default value

        const FloatT visibility_probability_0 = 0.5f;
        const FloatT visibility_probability_1 = visibilityProbability(color_similarity_ncc);

        const FloatT message_prob_0 = prev_backward_message * visibility_probability_1 *    change_prob
                + (1 - prev_backward_message) * visibility_probability_0 * no_change_prob;
        const FloatT message_prob_1 = prev_backward_message * visibility_probability_1 * no_change_prob
                + (1 - prev_backward_message) * visibility_probability_0 *    change_prob;

        return message_prob_1 / (message_prob_0 + message_prob_1);
    }

    FloatT visibilityProbability(const FloatT color_similarity)
    {
        return std::exp(- (1 - color_similarity) * (1 - color_similarity) * inv_sigma_sq) * visibility_probability_normalization;
    }

private:
    static FloatT VisibilityProbabilityNormalization(const FloatT sigma)
    {
        // http://www.wolframalpha.com/input/?i=integrate+exp(-(1-x)%5E2+%2F+(2*s%5E2))+dx+from+-1+to+1
        return std::sqrt(M_PI / 2) * sigma * std::erf(std::sqrt(2.0) / sigma);
    }

public:
    const FloatT inv_sigma_sq;
    const FloatT visibility_probability_normalization;

    const FloatT no_change_prob, change_prob;
};
