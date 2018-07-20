#include "pwmvs.hpp"

#include "geometric_priors.hpp"
#include "random_generators.hpp"


#include "openMVG/image/image_converter.hpp"
#include "openMVG/image/image_io.hpp"
#include "image_io.hpp"

void DebugMinMaxDepth(const std::shared_ptr<RefView> &ref)
{
    std::cout << "ref.depth.min: " << ref->depth.minCoeff() << " max: " << ref->depth.maxCoeff() << std::endl;
}

PWMVS::PWMVS(std::shared_ptr<RefView> &ref, const std::vector<std::shared_ptr<SrcView> > &srcs, const PWMVS::Options &options)
    : ref(ref)
    , srcs(srcs)
    , options(options)
    , photometric_consistency_calculator(*ref, options.window_size, options.color_dispersion, options.spatial_dispersion)
    , geometric_consistency_calculator(*ref, options.geometric_consistency_max_error)
    , message_calculator(options.color_similarity_sigma)
    , rng(deg2rad(80))
{
    assert(srcs.size() > 0);

    ref->depth.resize(ref->width, ref->height);
    ref->normal.resize(ref->width, ref->height, true, Normal(0, 0, 0));

    //DebugDepth(ref, std::string("depth_initial") + std::string(".jpeg"));
    //DebugNormal(ref, std::string("normal_initial") + std::string(".jpeg"));

    initializeRefView();
    createMatrices();

    //DebugDepth(ref, std::string("depth_random") + std::string(".jpeg"));
    //DebugNormal(ref, std::string("normal_random") + std::string(".jpeg"));
}

void PWMVS::run()
{
    computeInitialColorSimilarities();

    DirectionTB top_to_bottom(ref->height, ref->width);
    DirectionRL right_to_left(ref->height, ref->width);
    DirectionBT bottom_to_top(ref->height, ref->width);
    DirectionLR left_to_right(ref->height, ref->width);

    Direction * directions[4] = { &top_to_bottom, &right_to_left, &bottom_to_top, &left_to_right };

    for (int iter = 0; iter < options.num_iterations; iter++)
    {
        for (int dir = 0; dir < 4; dir++)
        {
            //const FloatT perturbation = 0.5 / std::pow(2.0, iter + dir / 4.0);
            const FloatT perturbation = options.perturbation / std::pow(2.0, iter + dir / 4.0);
            const FloatT state_transition_probability = StateTransitionProbability(iter * 4 + dir, options.num_iterations * 4);

            std::cout << "iter:" << iter << " dir:" << dir << std::endl;
            std::cout << "perturbation:" << perturbation << " state_transition_probability:" << state_transition_probability << std::endl;

            sweep(perturbation, state_transition_probability, *directions[dir]);

            // DEBUG
            //DebugDepth(ref, std::string("depth_") + std::to_string(iter) + std::string("_") + std::to_string(dir) + std::string(".jpeg"), options.min_depth, options.max_depth);
            //DebugNormal(ref, std::string("normal_") + std::to_string(iter) + std::string("_") + std::to_string(dir) + std::string(".jpeg"));
            //DebugNormalIncident(ref, std::string("incident_") + std::to_string(iter) + std::string("_") + std::to_string(dir) + std::string(".jpeg"));

            prev_selection_probabilities.swap(selection_probabilities);
        }
    }

    filter();

    // DEBUG
    //DebugDepth(ref, std::string("depth_filtered") + std::string(".jpeg"), options.min_depth, options.max_depth);
    //DebugNormal(ref, std::string("normal_filtered") + std::string(".jpeg"));

    for (auto iter : hypothesis_count)
    {
        std::cout << "hypothesis: " << iter.first << " -> " << iter.second << std::endl;
    }
}

void PWMVS::sweep(const FloatT perturbation, const FloatT state_transition_probability, const PWMVS::Direction &dir)
{
    Image<FloatT> forward_messages(dir.dim1(), srcs.size(), true, 0.5);

    // Initialize seeds
    std::vector<std::mt19937::result_type> seeds(dir.dim1());
    for (int dim1 = 0; dim1 < dir.dim1(); dim1++)
        seeds[dim1] = urng();

    #pragma omp parallel for schedule(dynamic)
    for (int dim1 = options.window_size; dim1 < (dir.dim1() - options.window_size); dim1++)
    {
        decltype(urng) local_urng(seeds[dim1]);

        // compute backwards message
        for (int src_id = 0; src_id < srcs.size(); src_id++)
        {
            FloatT previous_backward_message = 0.5;

            for (int dim2 = dir.dim2() - 1 - options.window_size; dim2 >= options.window_size ; dim2--)
            {
                int row, col;
                dir(dim1, dim2, row, col);

                backward_messages[src_id](row, col) = previous_backward_message = message_calculator.backwardMessage(color_similarities[src_id](row, col), previous_backward_message);
            }
        }

        // compute selection probabilites
        for (int dim2 = options.window_size; dim2 < (dir.dim2() - options.window_size); dim2++)
        {
            int row, col;
            dir(dim1, dim2, row, col);

            int prev_row, prev_col;
            dir(dim1, std::max(options.window_size, dim2 - 1), prev_row, prev_col);

            const Vector2i curr(col, row);
            const Vector2i prev(prev_col, prev_row);

            FloatT current_depth = ref->depth(curr);
            Normal current_normal = ref->normal(curr);

            FloatT prev_depth = ref->depth(prev);
            Normal prev_normal = ref->normal(prev);

            if (std::abs(current_normal.norm() - 1) > 0.1)
            {
                std::cerr << "error invalid normal at " << "row:" << row << " col:" << col << " = " << current_normal.transpose() << std::endl;
            }

            std::vector<FloatT> probabilities(srcs.size());
            for (int src_id = 0; src_id < srcs.size(); src_id++)
            {
                const FloatT color_similarity = color_similarities[src_id](row, col); // from previous sweep

                // Compute forward message
                FloatT forward_message = forward_messages(src_id, dim1) = message_calculator.forwardMessage(color_similarity, forward_messages(src_id, dim1));
                FloatT selection_probability = SelectionProbability(forward_message, backward_messages[src_id](curr), prev_selection_probabilities[src_id](curr), state_transition_probability);
                FloatT geometric_prior_probability = PriorProbability(*ref, *srcs[src_id], current_normal, current_depth, curr, options.window_size, options.incident_prior_sigma, options.min_triangulation_prior);

                assert(!std::isnan(forward_message));
                assert(!std::isnan(selection_probability));
                assert(!std::isnan(geometric_prior_probability));

                assert(forward_message >= 0);
                assert(selection_probability >= 0);
                assert(geometric_prior_probability >= 0);

                if (std::isnan(geometric_prior_probability) || geometric_prior_probability < 0)
                    std::cerr << "geometric_prior_probability: " << geometric_prior_probability << " row:" << row << " col:" << col << " src_id:" << src_id << std::endl;

                if (std::isnan(selection_probability) || selection_probability < 0)
                    std::cerr << "selection_probability: " << selection_probability << " row:" << row << " col:" << col << " src_id:" << src_id << std::endl;

                probabilities[src_id] = selection_probability * geometric_prior_probability;
            }

            const Normal ray = ref->ray(convertToFloat(curr));

            FloatT prp_depth  = PropagateDepth(*ref, convertToFloat(prev), convertToFloat(curr), prev_normal, prev_depth);
            FloatT prt_depth  = rdg.perturb(current_depth, perturbation, local_urng);
            Normal prt_normal = rng.perturb(current_normal, ray, perturbation, local_urng);
            FloatT rnd_depth  = rdg(options.min_depth, options.max_depth, local_urng);
            Normal rnd_normal = rng(-ray, local_urng);

            Normal avg_normal = (current_normal + prev_normal).normalized();
            FloatT avg_depth  = (current_depth + prev_depth) / 2;

            // TODO: why is this needed? [normals on the edge are skewed]
            prp_depth = std::max(std::min(prp_depth, options.max_depth), options.min_depth);
            prt_depth = std::max(std::min(prt_depth, options.max_depth), options.min_depth);
            rnd_depth = std::max(std::min(rnd_depth, options.max_depth), options.min_depth);

            Hypothesis hypotheses[] = {
                { current_depth, current_normal },  /* OK 1*/
                { prp_depth, prev_normal },         /* OK 2*/
//                { rnd_depth, current_normal },      /* PAPER */
                { current_depth, rnd_normal },      /* PAPER */
//                { rnd_depth, rnd_normal },          /* PAPER */
                { avg_depth, avg_normal },          /* MY */
                { prt_depth, current_normal },      /* OK 5*/
                { current_depth, prt_normal },      /* OK 4*/
                { prt_depth, prt_normal }           /* OK 3*/
            };

            Image<FloatT> hypothesis_photoconsistency_cache;
            const int best_hypothesis_idx = bestHypothesis(hypotheses, probabilities, curr, local_urng, hypothesis_photoconsistency_cache);
            const Hypothesis &best_hypothesis = hypotheses[best_hypothesis_idx];

#pragma omp critical
            hypothesis_count[best_hypothesis_idx]++;

            ref->depth(curr) = best_hypothesis.depth;
            ref->normal(curr) = best_hypothesis.normal;

            // recompute color similarities, forward messages and selection probability
            for (int src_id = 0; src_id < srcs.size(); src_id++)
            {
                FloatT color_similarity;
                if (best_hypothesis_idx == 0)
                {
                    color_similarity = color_similarities[src_id](curr);
                }
                else
                {
                    color_similarity = hypothesis_photoconsistency_cache(src_id, best_hypothesis_idx);
                    if (std::isnan(color_similarity))
                    {
                        color_similarity = photometric_consistency_calculator(*srcs[src_id], curr);
                    }
                    color_similarities[src_id](curr) = color_similarity;
                }

                const FloatT forward_message = forward_messages(src_id, dim1) = message_calculator.forwardMessage(color_similarity, forward_messages(src_id, dim1));
                const FloatT selection_probability = selection_probabilities[src_id](curr) = SelectionProbability(forward_message, backward_messages[src_id](curr), prev_selection_probabilities[src_id](curr), state_transition_probability);

                assert(!std::isnan(forward_message));
                assert(!std::isnan(selection_probability));

                assert(forward_message > 0);
                assert(selection_probability > 0);
            }
        }
    }
}

void PWMVS::initializeRefView()
{
    std::cout << __FUNCTION__ << " - begin" << std::endl;    // TODO: do not initialize if this is not a first run
    // TODO: use member random generators

    DebugMinMaxDepth(ref);

    ViewRandomDepthInitialization(*ref, options.min_depth, options.max_depth, rdg, urng);
    ViewRandomNormalInitialization(*ref, rng, urng);

    DebugMinMaxDepth(ref);

    std::cout << __FUNCTION__ << " - end" << std::endl;
}

void PWMVS::createMatrices()
{
    for (int src_id = 0; src_id < srcs.size(); src_id++)
    {
        color_similarities.emplace_back(Image<FloatT>(ref->width, ref->height, true, 0));

        backward_messages.emplace_back(Image<FloatT>(ref->width, ref->height, false));

        prev_selection_probabilities.emplace_back(Image<FloatT>(ref->width, ref->height, true, 0.5));
        selection_probabilities.emplace_back(Image<FloatT>(ref->width, ref->height, true, 0.5));
    }
}

void PWMVS::computeInitialColorSimilarities()
{
    std::cout << __FUNCTION__ << " - begin" << std::endl;
    for (int src_id = 0; src_id < srcs.size(); src_id++)
    {
        Image<FloatT> &color_similarity = color_similarities[src_id];

        for (int row = options.window_size; row < (ref->height - options.window_size); row++)
        {
#pragma omp parallel for
            for (int col = options.window_size; col < (ref->width - options.window_size); col++)
            {
                if (ref->depth(row, col) != 0)
                {
                    color_similarity(row, col) = photometric_consistency_calculator(*srcs[src_id], Vector2i(col, row));
                }
            }
        }

        // DEBUG COLOR SIMILARITIES
        //DebugImage(color_similarities[src_id], std::string("color_similarities_") + std::to_string(src_id) + std::string(".jpeg"));
        //std::cout << "initial min/max:" << color_similarities[src_id].minCoeff() << " " << color_similarities[src_id].maxCoeff() << std::endl;
    }
}

void PWMVS::filter()
{
    std::cout << __FUNCTION__ << " - begin" << std::endl;

    const FloatT min_visibility_probability = message_calculator.visibilityProbability(options.filter_min_color_similarity);

    for (int row = 0; row < ref->height; row++)
    {
//#pragma omp parallel for
        for (int col = 0; col < ref->width; col++)
        {
            const Vector2i x(col, row);
            const Vector3 X = ref->unproject(x);
            const Normal &n = ref->normal(x);

            int num_consistent = 0;
            for (int src_id = 0; src_id < srcs.size(); src_id++)
            {
                const FloatT triangulation_angle = TriangulationAngle(*srcs[src_id], X);
                if (triangulation_angle <= options.filter_min_triangulation_angle)
                    continue;

                const FloatT incident_angle = IncidentAngle(*srcs[src_id], X, n);
                if (incident_angle >= deg2rad(90))
                    continue;

                // TODO: switch to second version below. This is temporary until `selection_probabilities` aren't correct.
                //if (options.filter_photometric_consistency && (color_similarities[src_id](x) < options.filter_min_color_similarity))
                if (options.filter_photometric_consistency && (selection_probabilities[src_id](x) < min_visibility_probability))
                    continue;

                if (options.filter_geometric_consistency && geometric_consistency_calculator(*srcs[src_id], convertToFloat(x), X) >= options.filter_geometric_consistency_max_error)
                    continue;

                num_consistent++;
            }

            if (num_consistent < options.filter_min_num_consistent)
            {
                ref->depth(x) = 0;
                ref->normal(x) = Normal::Zero();
            }
        }
    }

    std::cout << __FUNCTION__ << " - end" << std::endl;

}

FloatT PWMVS::evaluateHypothesis(const SrcView *src, const Vector2i &x, const PWMVS::Hypothesis &hypothesis)
{
    return photometric_consistency_calculator(*src, x, hypothesis.normal, hypothesis.depth);
}
