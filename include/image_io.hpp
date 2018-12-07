#pragma once

#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>

#include "types.hpp"
#include "geometry.hpp"

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/image/image_converter.hpp"
#include "openMVG/image/image_filtering.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/image/pixel_types.hpp"
#include "openMVG/image/sample.hpp"

#include "stlplus3/filesystemSimplified/file_system.hpp"

namespace openMVG
{
namespace image
{
    template<>
    inline void Convert( const openMVG::image::RGBColor& valin, FloatT& out )
    {
        out = (0.3 * valin.r() + 0.59 * valin.g() + 0.11 * valin.b()) / 255.0;
    }
    template<>
    inline void Convert( const Normal& valin, openMVG::image::RGBColor& out )
    {
        for (int i = 0; i < 3; i++)
            out(i) = valin(i) * 256.0;
    }
    template<>
    inline void Convert( const openMVG::image::RGBColor& valin, Normal& out )
    {
        for (int i = 0; i < 3; i++)
            out(i) = valin(i) / 255.0;
    }
    template<>
    inline void Convert( const FloatT& valin, unsigned char& out )
    {
        out = std::min(static_cast<unsigned char>(valin * 256.0), static_cast<unsigned char>(255));
    }
}
}
static void GetRed(Image<openMVG::image::RGBColor> &image, Image<unsigned char> &red)
{
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> temp = Eigen::Map<Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic,3>>(reinterpret_cast<unsigned char *>(image.data()) + 0, image.Height(), image.Width(), Eigen::Stride<Eigen::Dynamic,3>(3 * image.Width(), 3));
    red.swap(temp);
}
static void GetGreen(Image<openMVG::image::RGBColor> &image, Image<unsigned char> &green)
{
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> temp = Eigen::Map<Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic,3>>(reinterpret_cast<unsigned char *>(image.data()) + 1, image.Height(), image.Width(), Eigen::Stride<Eigen::Dynamic,3>(3 * image.Width(), 3));
    green.swap(temp);
}
static void GetBlue(Image<openMVG::image::RGBColor> &image, Image<unsigned char> &blue)
{
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> temp = Eigen::Map<Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic,3>>(reinterpret_cast<unsigned char *>(image.data()) + 2, image.Height(), image.Width(), Eigen::Stride<Eigen::Dynamic,3>(3 * image.Width(), 3));
    blue.swap(temp);
}
static void SetRed(Image<openMVG::image::RGBColor> &image, Image<unsigned char> &red)
{
    Eigen::Map<Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic,3>>(reinterpret_cast<unsigned char *>(image.data()) + 0, image.Height(), image.Width(), Eigen::Stride<Eigen::Dynamic,3>(3 * image.Width(), 3)) = red;
}
static void SetGreen(Image<openMVG::image::RGBColor> &image, Image<unsigned char> &green)
{
    Eigen::Map<Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic,3>>(reinterpret_cast<unsigned char *>(image.data()) + 1, image.Height(), image.Width(), Eigen::Stride<Eigen::Dynamic,3>(3 * image.Width(), 3)) = green;
}
static void SetBlue(Image<openMVG::image::RGBColor> &image, Image<unsigned char> &blue)
{
    Eigen::Map<Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic,3>>(reinterpret_cast<unsigned char *>(image.data()) + 2, image.Height(), image.Width(), Eigen::Stride<Eigen::Dynamic,3>(3 * image.Width(), 3)) = blue;
}

template <typename Image>
static void Undistort(
    const Image& imageIn,
    const openMVG::cameras::IntrinsicBase *cam,
    Image &image_ud,
    typename Image::Tpixel fillcolor = typename Image::Tpixel(0))
{
    // no distortion, perform a direct copy
    if (!cam->have_disto())
    {
        image_ud = imageIn;
        return;
    }

    // There is distortion
    image_ud.resize(imageIn.Width(), imageIn.Height(), true, fillcolor);
    const openMVG::image::Sampler2d<openMVG::image::SamplerCubic> sampler;

#pragma omp parallel for
    for (int row = 0; row < imageIn.Height(); ++row)
    {
        for (int col = 0; col < imageIn.Width(); ++col)
        {
            const Vector2i undisto_pix(col, row);
            // compute coordinates with distortion
            const Vector2 disto_pix = cam->get_d_pixel(convertToFloat(undisto_pix).cast<double>()).cast<FloatT>();
            // pick pixel if it is in the image domain
            if (imageIn.Contains(disto_pix.y(), disto_pix.x()))
            {
                image_ud(row, col) = sampler(imageIn, convertToOpenMVG(disto_pix.y()), convertToOpenMVG(disto_pix.x()));
            }
        }
    }
}

template <typename Image>
static void TransformImage(
    const Image& image_in,
    const openMVG::cameras::IntrinsicBase *cam_in,
    Image &image_out,
    const openMVG::cameras::IntrinsicBase *cam_out,
    typename Image::Tpixel fillcolor = typename Image::Tpixel(0))
{
    // There is distortion
    image_out.resize(cam_out->w(), cam_out->h(), true, fillcolor);
    const openMVG::image::Sampler2d<openMVG::image::SamplerCubic> sampler;

    Image temp;
    const Image *src = &image_in;

    if (cam_out->w() != cam_in->w() || cam_out->h() != cam_in->h())
    {
        // TODO: Use better method to compute proper amount of blur needed.
        const FloatT inv_scale_x = static_cast<FloatT>(cam_in->w()) / static_cast<FloatT>(cam_out->w());
        const FloatT inv_scale_y = static_cast<FloatT>(cam_in->h()) / static_cast<FloatT>(cam_out->h());
        const FloatT sigma_x = 1.05 * std::sqrt(std::max(static_cast<FloatT>(0), inv_scale_x - 1));
        const FloatT sigma_y = 1.05 * std::sqrt(std::max(static_cast<FloatT>(0), inv_scale_y - 1));

        openMVG::image::ImageGaussianFilter(image_in, (sigma_x + sigma_y) / 2, temp);
        src = &temp;
    }

    #pragma omp parallel for
    for (int row = 0; row < image_out.Height(); ++row)
    {
        for (int col = 0; col < image_out.Width(); ++col)
        {
            const Vector2 x = convertToFloat(Vector2i(col, row));
            const Vector2 x_ = cam_in->cam2ima(cam_in->add_disto(cam_out->remove_disto(cam_out->ima2cam(x.cast<double>())))).cast<FloatT>();

            // pick pixel if it is in the image domain
            if (image_in.Contains(x_.y(), x_.x()))
            {
                image_out(row, col) = sampler(*src, convertToOpenMVG(x_.y()), convertToOpenMVG(x_.x()));
            }
        }
    }
}

template <>
void TransformImage(
        const Image<openMVG::image::RGBColor>& image_in,
        const openMVG::cameras::IntrinsicBase *cam_in,
        Image<openMVG::image::RGBColor> &image_out,
        const openMVG::cameras::IntrinsicBase *cam_out,
        openMVG::image::RGBColor fillcolor)
{
    Image<openMVG::image::RGBColor> temp = image_in;

    Image<unsigned char> red_in, green_in, blue_in;
    Image<unsigned char> red_out, green_out, blue_out;
    GetRed(temp, red_in);
    GetGreen(temp, green_in);
    GetBlue(temp, blue_in);

    TransformImage(red_in, cam_in, red_out, cam_out, fillcolor);
    TransformImage(green_in, cam_in, green_out, cam_out, fillcolor);
    TransformImage(blue_in, cam_in, blue_out, cam_out, fillcolor);

    image_out.resizeLike(red_out);

    SetRed(image_out, red_out);
    SetGreen(image_out, green_out);
    SetBlue(image_out, blue_out);
}

static void NormalizeDepth(const Image<FloatT> &depth, const FloatT min_val, const FloatT max_val, Image<FloatT> &out)
{
    if (min_val == max_val)
    {
        out = depth;
        return;
    }

    out.resizeLike(depth);
    for (int row = 0; row < depth.Height(); row++)
        for (int col = 0; col < depth.Width(); col++)
            out(row, col) = std::max(static_cast<FloatT>(0), std::min(static_cast<FloatT>(1), (depth(row, col) - min_val) / (max_val - min_val)));
}

static void NormalizeNormals(const Image<Normal> &depth, Image<Normal> &out)
{
    out.resizeLike(depth);
    for (int row = 0; row < depth.Height(); row++)
        for (int col = 0; col < depth.Width(); col++)
        {
            Normal n = depth(row, col);
            n.array() += 1;
            n.array() /= 2;
            out(row, col) = n;
        }
}

template <class T>
static void NormalizeDepth(const Image<T> &depth, Image<T> &out)
{
    const FloatT min_val = depth.minCoeff();
    const FloatT max_val = depth.maxCoeff();

    NormalizeDepth(depth, min_val, max_val, out);
}

template <class T>
static bool Load(const std::string &filename, Image<T> &image)
{
    {
        Image<openMVG::image::RGBColor> tmp;
        if (openMVG::image::ReadImage(filename.c_str(), &tmp))
        {
            openMVG::image::ConvertPixelType(tmp, &image);
            return true;
        }
    }
    {
        Image<unsigned char> tmp;
        if (openMVG::image::ReadImage(filename.c_str(), &tmp))
        {
            openMVG::image::ConvertPixelType(tmp, &image);
            return true;
        }
    }

    std::cerr << "failed to load image: " << filename << std::endl;
    return false;
}

static bool StringCompare(const std::string &a, const std::string &b)
{
    if (a.size() != b.size())
        return false;

    for (auto a_iter = a.begin(), b_iter = b.begin(); a_iter != a.end() && b_iter != b.end(); ++a_iter, ++b_iter)
    {
        if (::tolower(*a_iter) != ::tolower(*b_iter))
            return false;
    }

    return true;
}

template <class T>
static bool Write(const std::string &filename, const Image<T> &image)
{
    const std::string ext = stlplus::extension_part(filename);
    if (StringCompare(ext, "jpg") || StringCompare(ext, "jpeg"))
    {
        return openMVG::image::WriteJpg(filename.c_str(), image, 100);
    }
    else if (!StringCompare(ext, "bin"))
    {
        return openMVG::image::WriteImage(filename.c_str(), image);
    }

    auto begin = reinterpret_cast<const char *>(image.data());
    auto end = reinterpret_cast<const char *>(image.data() + image.Width() * image.Height());

    std::ofstream out(filename.c_str(), std::ios::binary);
    if (!out.is_open())
    {
        std::cout << "failed to open file: " << filename << std::endl;
        return false;
    }

    unsigned char separator = ';';
    size_t width = image.Width(), height = image.Height();

    out << width << separator << height << separator << sizeof(T) << separator;

    std::copy(begin, end, std::ostreambuf_iterator<char>(out));
    out.flush();

    return !out.fail();
}

template <class T>
static bool Read(const std::string &filename, Image<T> &image)
{
    if (stlplus::extension_part(filename) != std::string("bin"))
    {
        return Load(filename, image);
    }

    std::ifstream in(filename.c_str(), std::ios::binary);
    if (!in.is_open())
    {
        std::cerr << "failed to open file: " << filename << std::endl;
        return false;
    }

    size_t width, height, depth;
    unsigned char separator;
    in >> width >> separator >> height >> separator >> depth >> separator;

    if (depth != sizeof(T))
    {
        std::cerr << "failed to read file (invalid depth): " << filename << std::endl;
        return false;
    }

    image.resize(width, height);

    in.read(reinterpret_cast<char *>(image.data()), image.Width() * image.Height() * sizeof(T));

    if (in.fail())
    {
        std::cerr << "failed to read file: " << filename << std::endl;
        return false;
    }

    return true;
}

static void DebugImage(const Image<FloatT> &image, const std::string &filename, FloatT min_depth, FloatT max_depth)
{
    Image<FloatT> depth_clone;
    Image<unsigned char> depth_unsigned;
    NormalizeDepth(image, min_depth, max_depth, depth_clone);
    openMVG::image::ConvertPixelType(depth_clone, &depth_unsigned);
    if (!openMVG::image::WriteImage(filename.c_str(), depth_unsigned))
        std::cerr << "error writing:" << filename << std::endl;
}

static void DebugImage(Image<FloatT> &image, const std::string &filename)
{
    DebugImage(image, filename, image.minCoeff(), image.maxCoeff());
}

static void DebugDepth(RefView *ref, const std::string &filename, FloatT min_depth, FloatT max_depth)
{
    DebugImage(ref->depth, filename, min_depth, max_depth);
}
static void DebugDepth(RefView *ref, const std::string &filename)
{
    DebugImage(ref->depth, filename);
}

static void DebugNormal(const Image<Normal> &normal, const std::string &filename)
{
    Image<Normal> normal_clone;
    Image<openMVG::image::RGBColor> normal_unsigned;
    NormalizeNormals(normal, normal_clone);
    openMVG::image::ConvertPixelType(normal_clone, &normal_unsigned);
    openMVG::image::WriteImage(filename.c_str(), normal_unsigned);
}
static void DebugNormal(RefView *ref, const std::string &filename)
{
    DebugNormal(ref->normal, filename);
}
static void DebugNormalIncident(RefView *ref, const std::string &filename)
{
    Image<FloatT> incident(ref->width, ref->height);

    for (int row = 0; row < ref->height; ++row)
    {
        for (int col = 0; col < ref->width; ++col)
        {
            const Vector2i x(col, row);
            incident(x) = std::min(deg2rad(90), angleBetweenVectors(-ref->ray(convertToFloat(x)), ref->normal(x))) / deg2rad(90);
        }
    }

    Image<unsigned char> incident_unsigned;
    openMVG::image::ConvertPixelType(incident, &incident_unsigned);
    openMVG::image::WriteImage(filename.c_str(), incident_unsigned);
}

static void ExportPoints(const View &v, const std::string &filename, const Matrix3 &R, const Vector3 &t)
{
    std::ofstream out(filename, std::ios::out | std::ios::trunc);

    long points = 0;
    for (int row = 0; row < v.height; row++)
        for (int col = 0; col < v.width; col++)
            if (v.depth(row, col) != 0)
                ++points;

    out << "ply" << "\n"
        << "format ascii 1.0" << "\n"
        << "comment calculated with PWMVS" << "\n"
        << "element vertex " << points << "\n"
        << "property float x" << "\n"
        << "property float y" << "\n"
        << "property float z" << "\n"
        << "property char red" << "\n"
        << "property char green" << "\n"
        << "property char blue" << "\n"
        << "property float nx" << "\n"
        << "property float ny" << "\n"
        << "property float nz" << "\n"
        << "end_header" << "\n";

    for (int row = 0; row < v.height; row++)
    {
        for (int col = 0; col < v.width; col++)
        {
            if (v.depth(row, col) != 0)
            {
                Vector2i x(col, row);
                Vector3 X = R * v.unproject(x) + t;
                int intensity = static_cast<char>(v.image(x) * 256);
                Normal n = R * v.normal(x);
                out << X.x() << " " << X.y() << " " << X.z() << " " << intensity << " " << intensity << " " << intensity << " " << n.x() << " " << n.y() << " " << n.z() << "\n";
            }
        }
    }
}

static void ExportPoints(const View &v, const std::string &filename)
{
    ExportPoints(v, filename, Matrix3::Identity(), Vector3::Zero());
}

static void ExportPoints(const std::vector<Vector3> &X, const std::vector<Normal> &n, const std::vector<openMVG::image::RGBColor> &c, const std::string &filename)
{
    std::ofstream out(filename, std::ios::out | std::ios::trunc);

    int size = X.size();

    out << "ply" << "\n"
        << "format ascii 1.0" << "\n"
        << "comment calculated with PWMVS" << "\n"
        << "element vertex " << size << "\n"
        << "property float x" << "\n"
        << "property float y" << "\n"
        << "property float z" << "\n"
        << "property uchar red" << "\n"
        << "property uchar green" << "\n"
        << "property uchar blue" << "\n"
        << "property float nx" << "\n"
        << "property float ny" << "\n"
        << "property float nz" << "\n"
        << "end_header" << "\n";

    const int max_char = std::numeric_limits<signed char>::max();

    for (int i = 0; i < size; ++i)
        out << X[i].x() << " " << X[i].y() << " " << X[i].z() << " " << (unsigned int)c[i].r() << " " << (unsigned int)c[i].g() << " " << (unsigned int)c[i].b() << " " << n[i].x() << " " << n[i].y() << " " << n[i].z() << "\n";
}


