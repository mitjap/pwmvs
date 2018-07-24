#pragma once

#include "types.hpp"

static FloatT cosAngleBetweenNormals(const Normal &a, const Normal &b)
{
    const FloatT cos_angle = a.dot(b);
    return std::min(std::max(cos_angle, static_cast<FloatT>(-1.0)), static_cast<FloatT>(+1.0));
}

static FloatT angleBetweenNormals(const Normal &a, const Normal &b)
{
    return std::acos(cosAngleBetweenNormals(a, b));
}

static FloatT cosAngleBetweenVectors(const Vector3 &a, const Vector3 &b)
{
    const FloatT cos_angle = a.dot(b) / (a.norm() * b.norm());
    return std::min(std::max(cos_angle, static_cast<FloatT>(-1.0)), static_cast<FloatT>(+1.0));
}

static FloatT angleBetweenVectors(const Vector3 &a, const Vector3 &b)
{
    return std::acos(cosAngleBetweenVectors(a, b));
}

static constexpr FloatT deg2rad(const FloatT deg)
{
    return deg * M_PI / 180;
}

static constexpr FloatT rad2deg(const FloatT rad)
{
    return rad * 180 / M_PI;
}

static FloatT QuadrilateralArea(const Vector2 &A, const Vector2 &B, const Vector2 &C, const Vector2 &D)
{
    // Area of polygon split to two triangles, polygon MUST always be convex
    return std::abs(0.5 * (
        A.x() * B.y() - B.x() * A.y() +
        D.x() * A.y() - A.x() * D.y() +
        B.x() * C.y() - C.x() * B.y() +
        C.x() * D.y() - D.x() * C.y()
    ));
}

// TODO: maybe this should be moved to types.hpp
static FloatT PropagateDepth(const RefView &view, const Vector2 &from, const Vector2 &to, const Normal &from_n, const FloatT from_depth)
{
    assert(std::abs(from_n.norm() - 1) < 1e-5);

    const Vector3 X_from = view.unproject(from, from_depth);
    Eigen::Hyperplane<FloatT, 3> plane(from_n, X_from);
    Eigen::ParametrizedLine<FloatT, 3> ray(Vector3::Zero(), view.ray(to));
    Vector3 X_to = ray.intersectionPoint(plane);

    const FloatT to_depth = view.distance(X_to);

    if (to_depth > 0)
        return to_depth;

    // If angle between `from_n` and `ray` is near 90° propagated point might be behind camera so just return initial depth
    return from_depth;
}

static Vector3 PerpendicularTo(const Vector3 &a) {
    int max_index; a.maxCoeff(&max_index);
    Vector3 b = Vector3::Zero();
    b((max_index + 1) % a.size()) = -a(max_index);
    b(max_index) = a((max_index + 1) % a.size());

    return b;
}

static Vector3 PerpendicularTo(const Vector3 &a, const Vector3 &b) {
    return a.cross(b);
}
