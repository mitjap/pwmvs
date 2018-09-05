#pragma once

#include "types.hpp"
#include "random_generators.hpp"
#include "photometric_consistency.hpp"
#include "geometric_consistency.hpp"
#include "geometric_priors.hpp"

class PWMVS
{

public:
    // There are many options, use struct to pass them
    struct Options
    {
        // Random depth initialization
        FloatT min_depth = 0;
        FloatT max_depth = std::numeric_limits<FloatT>::max();

        // General
        int window_size = 5;
        int num_iterations = 5;
        FloatT perturbation = 1;
        int monte_carlo_samples = 15;

        bool geometric_consistency_term = false;
        FloatT geometric_consistency_weight = 0.1;
        FloatT geometric_consistency_max_error = 5;

        // NCC weight parameters
        FloatT spatial_dispersion = 3.0;
        FloatT color_dispersion = 0.3;

        // Geometric priors
        FloatT incident_prior_sigma = deg2rad(45);
        FloatT min_triangulation_prior = deg2rad(5);

        // Message propagation (forward, backward)
        FloatT color_similarity_sigma = 0.6;


        // Filtering
        bool filter_photometric_consistency = true;
        FloatT filter_min_color_similarity = 0.1;
        FloatT filter_min_triangulation_angle = deg2rad(3.0);
        int filter_min_num_consistent = 2;

        bool filter_geometric_consistency = false;
        FloatT filter_geometric_consistency_max_error = 1;

    };

private:
    struct Hypothesis {
        FloatT depth;
        Normal normal;
    };

    struct Direction
    {
        Direction(int rows, int cols)
            : rows(rows)
            , cols(cols)
        {
            // nothing to do
        }

        virtual int dim1() const = 0;
        virtual int dim2() const = 0;
        virtual void operator()(int dim1, int dim2, int &row, int &col) const = 0;

    protected:
        int rows, cols;
    };

    struct DirectionTB : Direction
    {
        DirectionTB(int rows, int cols) : Direction(rows, cols) {}

        virtual int dim1() const override { return cols; }
        virtual int dim2() const override { return rows; }
        virtual void operator()(int dim1, int dim2, int &row, int &col) const override {
            row = dim2;
            col = dim1;
        }
    };

    struct DirectionBT : Direction
    {
        DirectionBT(int rows, int cols) : Direction(rows, cols) {}

        virtual int dim1() const override { return cols; }
        virtual int dim2() const override { return rows; }
        virtual void operator()(int dim1, int dim2, int &row, int &col) const override {
            row = rows - dim2 - 1;
            col = dim1;
        }
    };

    struct DirectionLR : Direction
    {
        DirectionLR(int rows, int cols) : Direction(rows, cols) {}

        virtual int dim1() const override { return rows; }
        virtual int dim2() const override { return cols; }
        virtual void operator()(int dim1, int dim2, int &row, int &col) const override {
            row = dim1;
            col = dim2;
        }
    };

    struct DirectionRL : Direction
    {
        DirectionRL(int rows, int cols) : Direction(rows, cols) {}

        virtual int dim1() const override { return rows; }
        virtual int dim2() const override { return cols; }
        virtual void operator()(int dim1, int dim2, int &row, int &col) const override {
            row = dim1;
            col = cols - dim2 - 1;
        }
    };

public:
    PWMVS(std::shared_ptr<RefView> &ref, const std::vector<std::shared_ptr<SrcView>> &srcs, const Options &options);

    bool run();
protected:

    void sweep(const FloatT perturbation, const FloatT state_transition_probability, const Direction &dir);
    void initializeRefView();
    void createMatrices();
    void computeInitialColorSimilarities();

    void filter();

    template <size_t N, class _UniformRandomNumberGenerator>
    int bestHypothesis(const Hypothesis (&hypotheses)[N], const std::vector<FloatT> &probabilities, const Vector2i &x, _UniformRandomNumberGenerator &urng, Image<FloatT> &hypothesis_photoconsistency_cache) {
        static_assert(N > 0, "N must be positive");

        // First hypothesis MUST be current solution (for optimization purposes)
        assert(hypotheses[0].depth == ref->depth(x));
        assert(hypotheses[0].normal == ref->normal(x));

        hypothesis_photoconsistency_cache.resize(N, srcs.size(), true, std::numeric_limits<FloatT>::quiet_NaN());

        FloatT evaluations[N];
        for (int i = 0; i < N; i++)
            evaluations[i] = 0;

        Sampler sampler(probabilities);

        // sample and accumulate evaluations per hypothesis
        for (int i = 0; i < options.monte_carlo_samples; i++) {
            int src_id = sampler(urng);
            const std::shared_ptr<SrcView> &src = srcs[src_id];

            evaluations[0] += color_similarities[src_id](x);
            for (int h = 1; h < N; ++h)
            {
                FloatT photo_consistency = hypothesis_photoconsistency_cache(src_id, h);
                if (std::isnan(photo_consistency))
                {
                    photo_consistency = photometric_consistency_calculator(*src, x, hypotheses[h].normal, hypotheses[h].depth);
                    hypothesis_photoconsistency_cache(src_id, h) = photo_consistency;
                }

                evaluations[h] += photo_consistency;
            }

            if (options.geometric_consistency_term)
            {
                for (int h = 0; h < N; ++h)
                {
                    // TODO: also consider src's normal at projected pixel
                    evaluations[h] += options.geometric_consistency_weight * geometric_consistency_calculator(*src, convertToFloat(x), hypotheses[h].depth);
                }
            }
        }

        // find best hypothesis
        int best_idx = 0;
        FloatT best_evaluation = evaluations[best_idx];
        for (int i = 1; i < N; i++)
        {
            const FloatT evaluation = evaluations[i];
            if (evaluation > best_evaluation)
            {
                best_idx = i;
                best_evaluation = evaluation;
            }
        }

        return best_idx;
    }

    FloatT evaluateHypothesis(const SrcView *src, const Vector2i &x, const Hypothesis &hypothesis);


private:
    const Options options;

    const std::shared_ptr<RefView> ref;
    const std::vector<std::shared_ptr<SrcView>> srcs;

    std::vector<Image<FloatT>> color_similarities; // NCC

    std::vector<Image<FloatT>> backward_messages;

    std::vector<Image<FloatT>> selection_probabilities;
    std::vector<Image<FloatT>> prev_selection_probabilities;

    //Image<FloatT> sampling_probabilities;   // first dimension is `image_id`, second dimension is image column
    //Image<FloatT> forward_messages;         // first dimension is `image_id`, second dimension is image column

    RandomDepthGenerator rdg;
    RandomNormalGenerator rng;

    PhotometricConsistencyOptimized<> photometric_consistency_calculator;
    GeometricConsistency geometric_consistency_calculator;
    MessageCalculator message_calculator;

    std::mt19937 urng;

    std::map<int, int> hypothesis_count;
};
