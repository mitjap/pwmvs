#pragma once

#include "types.hpp"

#include <openMVG/image/sample.hpp>

static Matrix3 ComputeHomography(const Matrix3 &ref_K_inv, const Matrix3 &src_K, const Matrix3 &src_R, const Vector3 &src_t, const Normal &n, const Vector3 &X)
{
    const FloatT d = n.dot(X);
    return src_K * (src_R + src_t * (n.transpose() / d)) * ref_K_inv;
}

static Matrix3 ComputeHomography(const RefView &ref, const SrcView &src, const Normal &n, const FloatT depth, const Vector2 &x)
{
    return ComputeHomography(ref.K_inv, src.K, src.R, src.t, n, ref.unproject(x, depth));
}

struct PhotometricConsistency
{
public: // TODO: for testing purposes
    struct AbstractWindow {
        virtual FloatT operator()(int row, int col) const = 0;
    };


    struct WindowProduct : public AbstractWindow
    {
        WindowProduct(const AbstractWindow &a, const AbstractWindow &b)
            : AbstractWindow()
            , a(a)
            , b(b)
        {
            // nothing to do
        }

        virtual FloatT operator()(int window_row, int window_col) const override
        {
            return a(window_row, window_col) * b(window_row, window_col);
        }

    private:
        const AbstractWindow &a;
        const AbstractWindow &b;
    };

    struct Window : public AbstractWindow {
        Window(const Image<FloatT> &image, const Vector2i &ref_center)
            : image(image)
            , ref_center(ref_center)
        {
            // nothing to do
        }


        WindowProduct operator *(const Window &other) const
        {
            return WindowProduct(*this, other);
        }

    protected:
        const Image<FloatT> &image;
        const Vector2i ref_center;
    };

    struct RefWindow : public Window {
        RefWindow(const Image<FloatT>& ref_image, const Vector2i &ref_center)
            : Window(ref_image, ref_center)
        {
            // nothing to do
        }

        virtual FloatT operator()(int window_row, int window_col) const override
        {
            int row = ref_center.y() + window_row;
            int col = ref_center.x() + window_col;

            if (!image.Contains(row, col))
            {
                row = std::max(image.Height() - 1, std::max(0, row));
                col = std::max(image.Width()  - 1, std::max(0, col));
                return image(row, col);
            }

            //return image(row, col);

            return sampler(image, convertToOpenMVG(convertToFloat(row)), convertToOpenMVG(convertToFloat(col)));
        }


        openMVG::image::Sampler2d<openMVG::image::SamplerLinear> sampler;
    };

    struct SrcWindow : public Window {
        SrcWindow(const Image<FloatT>& src_image, const Vector2i &ref_center, const Matrix3 &H)
            : Window(src_image, ref_center)
            , H(H)
            , src_center(transfer(0, 0))
        {
            // nothing to do
        }

        virtual FloatT operator()(int window_row, int window_col) const override
        {
            Vector2 x = transfer(window_row, window_col);

            if (!image.Contains(x))
            {
                x(0) = std::min(FloatT(image.Height()), std::max(FloatT(0), x(0)));
                x(1) = std::min(FloatT(image.Width()), std::max(FloatT(0), x(1)));
            }

            //return image(convertToInt(x));

            return sampler(image, convertToOpenMVG(x.y()), convertToOpenMVG(x.x()));
        }

    private:
        Vector2 transfer(int window_row, int window_col) const
        {
            return (H * convertToFloat(Vector2i(window_col + ref_center.x(), window_row + ref_center.y())).homogeneous()).hnormalized();
        }

    private:
        const Matrix3 &H;
        const Vector2 src_center;


        openMVG::image::Sampler2d<openMVG::image::SamplerLinear> sampler;
    };

public:

    /// NOTE: this implementation is not optimized at all. It serves as reference implementation. Many calculations repeat many times and could be cached.
    PhotometricConsistency(const RefView &ref, const int window_size, const FloatT color_dispersion, FloatT spatial_dispersion)
        : window_size(window_size)
        , color_dispersion(1 / (2 * color_dispersion * color_dispersion))
        , spatial_dispersion(1 / (2 * spatial_dispersion * spatial_dispersion))
        , ref(ref)
    {
        assert(window_size > 0);
        assert(color_dispersion > 0);
        assert(spatial_dispersion > 0);
    }


public: // TODO: for testing purposes


    FloatT operator ()(const SrcView &src, const Vector2i &x)
    {
        return (*this)(src, x, ref.normal(x), ref.depth(x));
    }

    FloatT operator ()(const SrcView &src, const Vector2i &x, const Normal &normal, FloatT depth)
    {
        if (!src.isVisible(src.project(ref.unproject(convertToFloat(x), depth)))) return -1;

        const Matrix3 H = ComputeHomography(ref, src, normal, depth, convertToFloat(x));

        RefWindow ref_window(ref.image, x);
        SrcWindow src_window(src.image, x, H);

        FloatT cov_src_src = covariance(ref_window, src_window, src_window);
        if (cov_src_src < 1e-5) return -1;
        FloatT cov_ref_ref = covariance(ref_window, ref_window, ref_window);
        if (cov_src_src < 1e-5) return -1;

        FloatT cov_ref_src = covariance(ref_window, ref_window, src_window);

        FloatT ncc = cov_ref_src / std::sqrt(cov_ref_ref * cov_src_src);
        ncc = std::max(static_cast<FloatT>(-1), std::min(static_cast<FloatT>(+1), ncc));
        return ncc;
    }

    FloatT covariance(const RefWindow &ref, const Window &a, const Window &b)
    {
        return weightedAverage(ref, a * b) - weightedAverage(ref, a) * weightedAverage(ref, b);
    }

    FloatT weightedAverage(const RefWindow &ref_window, const AbstractWindow &window)
    {
        FloatT ref_center_color = ref_window(0, 0);

        FloatT weighted_color_sum = 0;
        FloatT weight_sum = 0;
        for (int row = -window_size; row <= +window_size; row++)
        {
            for (int col = -window_size; col <= +window_size; col++)
            {
                FloatT color = window(row, col);
                FloatT w = weight(Vector2(col, row), ref_window(row, col) - ref_center_color);

                weighted_color_sum += w * color;
                weight_sum += w;
            }
        }

        return weighted_color_sum / weight_sum;
    }

    FloatT weight(const Vector2 &spatial_diff, const FloatT color_diff)
    {
        return std::exp(- (color_diff * color_diff) * color_dispersion - spatial_diff.squaredNorm() * spatial_dispersion);
    }

private:
    int window_size;
    FloatT color_dispersion;
    FloatT spatial_dispersion;

    const RefView &ref;
};

template <class Sampler = openMVG::image::SamplerLinear>
struct PhotometricConsistencyOptimized
{
    PhotometricConsistencyOptimized(const RefView &ref, const int window_size, const FloatT color_dispersion, FloatT spatial_dispersion)
        : window_size(window_size)
        , color_dispersion(1 / (2 * color_dispersion * color_dispersion))
        , spatial_dispersion(1 / (2 * spatial_dispersion * spatial_dispersion))
        , ref(ref)
        , ref_weighted_color_sum(ref.image.Width(), ref.image.Height(), true, 0)
        , ref_weighted_color_sum_sq(ref.image.Width(), ref.image.Height(), true, 0)
    {
        assert(window_size > 0);
        assert(color_dispersion > 0);
        assert(spatial_dispersion > 0);

        precalculate();
    }

public:
    FloatT operator ()(const SrcView &src, const Vector2i &x)
    {
        return (*this)(src, x, ref.normal(x), ref.depth(x));
    }

    FloatT operator ()(const SrcView &src, const Vector2i &x, const Normal &normal, FloatT depth)
    {
        const Vector3 X = ref.unproject(convertToFloat(x), depth);

        if (!src.isVisible(src.project(X))) return -1;

        const FloatT center_color = Sample(ref.image, x);

        FloatT weights_sum = 0;
        FloatT weighted_src_color_sum = 0;
        FloatT weighted_src_ref_color_sum = 0;
        FloatT weighted_src_color_sum_sq = 0;

        const FloatT weighted_ref_color_sum = ref_weighted_color_sum(x);
        const FloatT weighted_ref_color_sum_sq = ref_weighted_color_sum_sq(x);

        const FloatT d = X.dot(normal);

        for (int window_row = -window_size; window_row <= window_size; window_row++)
        {
            for (int window_col = -window_size; window_col <= window_size; window_col++)
            {
                const Vector2i x_(window_col, window_row);

                const Normal ray = ref.ray(convertToFloat((x + x_).eval()));
                const Vector2 x_src = src.project((d / ray.dot(normal)) * ray);

                const FloatT ref_color = Sample(ref.image, (x + x_).eval());
                const FloatT src_color = Sample(src.image, x_src, false);

                const FloatT w = weight(window_row, window_col, center_color - ref_color);

                weighted_src_color_sum     += w * src_color;
                weighted_src_color_sum_sq  += w * src_color * src_color;
                weighted_src_ref_color_sum += w * src_color * ref_color;
                weights_sum += w;
            }
        }

        weighted_src_color_sum     /= weights_sum;
        weighted_src_color_sum_sq  /= weights_sum;
        weighted_src_ref_color_sum /= weights_sum;

        const FloatT ref_color_var = weighted_ref_color_sum_sq - (weighted_ref_color_sum * weighted_ref_color_sum);
        const FloatT src_color_var = weighted_src_color_sum_sq - (weighted_src_color_sum * weighted_src_color_sum);
        const FloatT ref_src_color_covar = weighted_src_ref_color_sum - (weighted_ref_color_sum * weighted_src_color_sum);

        if (ref_color_var < 1e-6) return -1;
        if (src_color_var < 1e-6) return -1;

        FloatT ncc = ref_src_color_covar / std::sqrt(ref_color_var * src_color_var);
        ncc = std::max(static_cast<FloatT>(-1), std::min(static_cast<FloatT>(+1), ncc));
        return ncc;
    }

private:
    void precalculate()
    {
        for (int row = window_size; row < (ref.image.Height() - window_size); row++)
        {
#pragma omp parallel for
            for (int col = window_size; col < (ref.image.Width() - window_size); col++)
            {
                const FloatT center_color = Sample(ref.image, Vector2i(col, row));

                FloatT weights_sum = 0;
                FloatT weighted_color_sum = 0;
                FloatT weighted_color_sum_sq = 0;

                for (int window_row = -window_size; window_row <= window_size; window_row++)
                {
                    for (int window_col = -window_size; window_col <= window_size; window_col++)
                    {
                        const FloatT color = Sample(ref.image, Vector2i(col + window_col, row + window_row));
                        const FloatT w = weight(window_row, window_col, center_color - color);

                        weighted_color_sum    += w * color;
                        weighted_color_sum_sq += w * color * color;
                        weights_sum += w;
                    }
                }

                ref_weighted_color_sum(row, col) = weighted_color_sum / weights_sum;
                ref_weighted_color_sum_sq(row, col) = weighted_color_sum_sq / weights_sum;
            }
        }
    }

    FloatT weight(int row_diff, int col_diff, const FloatT color_diff)
    {
        return std::exp(- (color_diff * color_diff) * color_dispersion - (row_diff * row_diff + col_diff * col_diff) * spatial_dispersion);
    }

    FloatT Sample(const Image<FloatT> &image, const Vector2 &x, bool clamp_to_edge) const
    {
        if (clamp_to_edge && !image.Contains(x))
        {
            Vector2i x_;
            x_(0) = std::min(image.Width() - 1, std::max(0, convertToInt(x(0))));
            x_(1) = std::min(image.Height() - 1, std::max(0, convertToInt(x(1))));

            return Sample(image, x_);
        }

        return sampler(image, convertToOpenMVG(x.y()), convertToOpenMVG(x.x()));
    }

    FloatT Sample(const Image<FloatT> &image, const Vector2i &x) const
    {
        return image(x);
    }

private:


private:
    int window_size;
    FloatT color_dispersion;
    FloatT spatial_dispersion;

    const RefView &ref;

    Image<FloatT> ref_weighted_color_sum;
    Image<FloatT> ref_weighted_color_sum_sq;

    openMVG::image::Sampler2d<Sampler> sampler;
};
