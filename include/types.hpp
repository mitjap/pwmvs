#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <random>

#include <openMVG/numeric/numeric.h>
#include <openMVG/image/image_container.hpp>
#include <openMVG/image/pixel_types.hpp>

#define DEBUG(var) std::cout << __FUNCTION__ << " " << __LINE__ << " " << #var" = " << var << std::endl

using FloatT = float;
using Normal = Eigen::Matrix<FloatT, 3, 1>;
using Matrix3 = Eigen::Matrix<FloatT, 3, 3>;
using Matrix34 = Eigen::Matrix<FloatT, 3, 4>;
using Vector3 = Eigen::Matrix<FloatT, 3, 1>;
using Vector2 = Eigen::Matrix<FloatT, 2, 1>;
using Vector2i = Eigen::Matrix<int, 2, 1>;


template <class T>
static int convertToInt(T x) {
    return std::round(x);

    //return static_cast<int>(x);
}

template <class T>
static float convertToFloat(T x) {
    return static_cast<float>(x);

    //return static_cast<float>(x) + 0.5f;
}

template <>
float convertToFloat(float x) {
    return x;
}

template <>
float convertToFloat(double x) {
    return static_cast<float>(x);
}


static Vector2 convertToFloat(const Vector2i &x) {
    Vector2 x_;
    x_(0) = convertToFloat(x(0));
    x_(1) = convertToFloat(x(1));
    return x_;
}

static Vector2i convertToInt(const Vector2 &x) {
    Vector2i x_;
    x_(0) = convertToInt(x(0));
    x_(1) = convertToInt(x(1));
    return x_;
}

template <class T>
T convertToOpenMVG(T val) {
    return val;

//    static_assert(std::is_floating_point<T>::value, "T not floating point type");
//    return val - static_cast<T>(0.5);
}
template <class T>
T convertFromOpenMVG(T val) {
    return val;

//    static_assert(std::is_floating_point<T>::value, "T not floating point type");
//    return val + static_cast<T>(0.5);
}

template<class T>
class Image : public openMVG::image::Image<T>
{
public:
    /// Full internal type
    using Base = openMVG::image::Image<T>;

    inline Image()
        : Base()
    {

    }

    Image( int width, int height, bool fInit = true, const T val = T() )
        : Base(width, height, fInit, val)
    {

    }

    using Base::Contains;

    bool Contains( const Vector2 &x ) const
    {
      return Contains(convertToInt(x));
    }

    bool Contains( const Vector2i &x ) const
    {
      return Base::Contains( x.y(), x.x() );
    }

    inline const T& operator()( int y, int x ) const
    {
      return Base::operator()( y, x );
    }

    inline T& operator()( int y, int x )
    {
      return Base::operator()( y, x );
    }

    inline const T& operator ()( const Vector2i &x) const
    {
      return Base::operator()( x.y(), x.x() );
    }

    inline T& operator ()( const Vector2i &x )
    {
      return Base::operator()( x.y(), x.x() );
    }

    inline const T& operator ()( const Vector2 &x) const
    {
      std::cerr << "WARNING: floating point acces to image" << std::endl;
      return Base::operator()( x.y(), x.x() );
    }

    inline T& operator ()( const Vector2 &x )
    {
      std::cerr << "WARNING: floating point acces to image" << std::endl;
      return Base::operator()( x.y(), x.x() );
    }
};


struct View
{
    View(int w, int h)
        : height(h)
        , width(w)
        , K(Matrix3::Identity())
        , K_inv(Matrix3::Identity())
    {
//        FloatT max_dim = std::max(width, height);
//        setIntrinsics(max_dim, max_dim/2, max_dim/2);
    }

//    void setIntrinsics(const FloatT focal, const FloatT ppx, const FloatT ppy)
//    {
//        assert(focal > 0);
//        assert(ppx > 0 && ppx < width);
//        assert(ppy > 0 && ppy < height);

//        K(0, 0) = K(1, 1) = focal;
//        K(0, 2) = ppx;
//        K(1, 2) = ppy;

//        K_inv = K.inverse();
//    }

    Normal ray(const Vector2 &x) const
    {
        return (K_inv * x.homogeneous()).normalized();
    }

    bool isVisible(const Vector2 &x) const {
        return x.x() >= 0 && x.x() < width
            && x.y() >= 0 && x.y() < height;
    }

    bool isVisible(const Vector2i &x) const {
        return x.x() >= 0 && x.x() < width
            && x.y() >= 0 && x.y() < height;
    }

    virtual Vector2 project(const Vector3 &X) const = 0;
    virtual Normal projectNormal(const Normal &n) const = 0;

    virtual Vector3 unproject(const Vector2 &x, const FloatT depth) const = 0;
    virtual Vector3 unproject(const Vector2i &x) const
    {
        assert(isVisible(x));
        return unproject(convertToFloat(x), depth(x));
    }

    virtual Normal unprojectNormal(const Normal &n) const = 0;
    Normal unprojectNormal(const Vector2i &x) const
    {
        return unprojectNormal(normal(x));
    }

    virtual FloatT distance(const Vector3 &X) const = 0;

    bool isValid(int row, int col) const
    {
        return depth(row, col) != 0;
    }

    bool isValid(const Vector2i &x) const
    {
        return depth(x) != 0;
    }

    Matrix3 K;
    Matrix3 K_inv;

    int width, height;

    Image<FloatT> image;
    Image<openMVG::image::RGBColor> color;
    Image<FloatT> depth;
    Image<Normal> normal;
};

struct SrcView : public View {
    SrcView(int width, int height)
        : SrcView(width, height, Matrix3::Identity(), Vector3::Zero())
    {
        // nothing to do
    }

    SrcView(int width, int height, const Matrix3 &R, const Vector3 &t)
        : View(width, height)
        , R(R)
        , t(t)
    {
        initialize();
    }

    void initialize()
    {
        C = -(R.transpose() * t);

        P.leftCols<3>() = K * R;
        P.rightCols<1>() = K * t;

        P_inv.leftCols<3>() = R.transpose() * K_inv;
        P_inv.rightCols<1>() = C;
    }


    virtual Vector2 project(const Vector3 &X) const override
    {
        return (P * X.homogeneous()).hnormalized();
    }

    using View::unproject;
    virtual Vector3 unproject(const Vector2 &x, const FloatT depth) const override
    {
        return P_inv * (depth * x.homogeneous()).homogeneous();
    }

    virtual Normal projectNormal(const Normal &n) const
    {
        return (R * n).normalized();
    }

    using View::unprojectNormal;
    virtual Normal unprojectNormal(const Normal &n) const
    {
        return (R.transpose() * n).normalized();
    }

    virtual FloatT distance(const Vector3 &X) const override
    {
        return P.bottomRows<1>() * X.homogeneous();
    }


    Matrix3 R;
    Vector3 t;

    Vector3 C;
    Matrix34 P;
    Matrix34 P_inv;
};

struct RefView : public View {
    RefView(int width, int height)
        : View(width, height)
    {
        // nothing to do
    }

    virtual Vector2 project(const Vector3 &X) const override
    {
        return (K * X).hnormalized();
    }

    using View::unproject;
    virtual Vector3 unproject(const Vector2 &x, const FloatT depth) const override
    {
        return depth * (K_inv * x.homogeneous());
    }

    virtual FloatT distance(const Vector3 &X) const override
    {
        return X.z();
    }

    virtual Normal projectNormal(const Normal &n) const
    {
        return n;
    }

    using View::unprojectNormal;
    virtual Normal unprojectNormal(const Normal &n) const
    {
        return n;
    }
};
