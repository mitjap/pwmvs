#pragma once

#include <random>

#include "types.hpp"
#include "geometry.hpp"

struct RandomDepthGenerator {
    RandomDepthGenerator()
    { }

    template<class UniformRandomNumberGenerator>
    FloatT operator()(const FloatT min_depth, const FloatT max_depth, UniformRandomNumberGenerator &urng) {
        assert(min_depth <= max_depth);

        std::uniform_real_distribution<FloatT> dist(min_depth, max_depth);
        return dist(urng);
    }

    template<class UniformRandomNumberGenerator>
    FloatT perturb(const FloatT depth, const FloatT perturbation, UniformRandomNumberGenerator &urng)
    {
        assert(0 <= perturbation && perturbation <= 1);

        FloatT min_depth = (1 - perturbation) * depth;
        FloatT max_depth = (1 + perturbation) * depth;

        return (*this)(min_depth, max_depth, urng);
    }
};


struct RandomNormalGenerator {
    RandomNormalGenerator(const FloatT max_angle)
        : dist(-1, +1)
        , dist_0_1(0, 1)
        , max_angle(max_angle)
        , cos_max_angle(cos(max_angle))
    { }

    // Random point on a sphere cap with maximum angle `max_angle` from `vec`
    template <class UniformRandomNumberGenerator>
    Normal operator()(const Vector3 &vec, UniformRandomNumberGenerator &urng) {
        return (*this)(vec, cos_max_angle, urng);
    }

    // Random point on a sphere cap with maximum angle `max_angle` from `vec`
    template <class UniformRandomNumberGenerator>
    Normal operator()(const Normal &a, const FloatT cos_max_angle, UniformRandomNumberGenerator &urng) {
        // Christian Blatter (https://math.stackexchange.com/users/1303/christian-blatter), Generate a random direction within a cone, URL (version: 2012-08-16): https://math.stackexchange.com/q/182936
        const Normal u = PerpendicularTo(a).normalized();
        const Normal v = PerpendicularTo(a, u).normalized();

        const FloatT phi = 2 * M_PI * dist_0_1(urng); // [0, 2 * M_PI)

        const FloatT z = dist_0_1(urng) * (1 - cos_max_angle) + cos_max_angle;
        const FloatT theta = std::acos(z); // [cos(max_angle), 1)

        const Normal n = std::sin(theta) * (std::cos(phi) * u + std::sin(phi) * v) + std::cos(theta) * a;

        assert(angleBetweenNormals(n, a) <= max_angle);
        return n;
    }

    template <class UniformRandomNumberGenerator>
    Normal perturb(const Normal &normal, const Normal &ray, FloatT perturbation, UniformRandomNumberGenerator &urng, int max_iter = 7)
    {
        const FloatT perturbation_angle = perturbation * (max_angle / 4);
        const Normal n = (*this)(normal, std::cos(perturbation_angle), urng);

        const FloatT cos_angle = cosAngleBetweenNormals(n, -ray);
        if (cos_angle >= cos_max_angle) {
            return n;
        } else {
            const FloatT angle = std::acos(cos_angle) - max_angle;
            const Normal axis = PerpendicularTo(n, -ray).normalized();

            Matrix3 R = Eigen::AngleAxis<FloatT>(2 * angle, axis).toRotationMatrix();
            return (R * n).normalized();
        }
    }

private:
    std::uniform_real_distribution<FloatT> dist;
    std::uniform_real_distribution<FloatT> dist_0_1;

public:
    const FloatT max_angle;
    const FloatT cos_max_angle;
};

template <class UniformRandomNumberGenerator>
static void ViewRandomDepthInitialization(View &v, const FloatT min_depth, const FloatT max_depth, RandomDepthGenerator &depth_gen, UniformRandomNumberGenerator &urng)
{
    // Initialize seeds
    std::vector<typename UniformRandomNumberGenerator::result_type> seeds(v.height);
    for (int row = 0; row < v.height; row++)
        seeds[row] = urng();

    #pragma omp parallel for
    for (int row = 0; row < v.height; row++)
    {
        UniformRandomNumberGenerator local_urng(seeds[row]);
        for (int col = 0; col < v.width; col++)
        {
            const FloatT d = depth_gen(min_depth, max_depth, local_urng);

            // some depth values might have already been set
            if (v.depth(row, col) == 0)
                v.depth(row, col) = d;
        }
    }
}

template <class UniformRandomNumberGenerator>
static void ViewRandomNormalInitialization(View &v, RandomNormalGenerator &normal_gen, UniformRandomNumberGenerator &urng)
{
    // Initialize seeds
    std::vector<typename UniformRandomNumberGenerator::result_type> seeds(v.height);
    for (int row = 0; row < v.height; row++)
        seeds[row] = urng();

    #pragma omp parallel for
    for (int row = 0; row < v.height; row++)
    {
        UniformRandomNumberGenerator local_urng(seeds[row]);
        for (int col = 0; col < v.width; col++)
        {
            const Vector2i x(col, row);

            Normal &n = v.normal(x);
            const Normal ray = v.ray(convertToFloat(x));

            if (std::abs(n.norm() - 1) >= std::numeric_limits<FloatT>::epsilon() // not a unit length vector
                || cosAngleBetweenNormals(n, -ray) < normal_gen.cos_max_angle)          // violates `max_angle`
            {
                n = normal_gen(-ray, local_urng);
            }

            assert(std::abs(n.norm() - 1) < 1e-3);
        }
    }
}

struct Sampler {
    Sampler(const std::vector<FloatT> &probabilities)
        : dist(probabilities.begin(), probabilities.end())
    { }

    template<class UniformRandomNumberGenerator>
    int operator()(UniformRandomNumberGenerator &urng) {
        return dist(urng);
    }

private:
    std::discrete_distribution<int> dist;
};

struct NaiveSampler {
    NaiveSampler(const std::vector<FloatT> &probabilities)
    {
        best_probability_idx = 0;
        for (int i = 1; i < probabilities.size(); i++)
        {
            if (probabilities[i] > probabilities[best_probability_idx])
                best_probability_idx = i;
        }
    }

    template<class UniformRandomNumberGenerator>
    int operator()(UniformRandomNumberGenerator &urng) {
        return best_probability_idx;
    }

private:
    int best_probability_idx;
};
