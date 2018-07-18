#pragma once

#include "types.hpp"

struct GeometricConsistency
{
    GeometricConsistency(const RefView &ref, const FloatT max_error)
        : ref(ref)
        , max_error(max_error)
    {
        // nothing to do
    }

    FloatT operator ()(const SrcView &src, const Vector2i &x_ref) const
    {
        return (*this)(src, convertToFloat(x_ref), ref.depth(x_ref));
    }

    FloatT operator ()(const SrcView &src, const Vector2 &x_ref, const FloatT depth) const
    {
        return (*this)(src, x_ref, ref.unproject(x_ref, depth));
    }
    FloatT operator ()(const SrcView &src, const Vector2 &x_ref, const Vector3 &X) const
    {
        const Vector2 x_src = src.project(X);

        if (!src.isVisible(x_src)) return max_error;

        // TODO: lineraly interpolate depth
        const Vector2 x_ref_ = ref.project(src.unproject(x_src, src.depth(convertToInt(x_src))));
        return std::min(max_error, (x_ref_ - x_ref).norm());
    }

    const RefView &ref;
    const FloatT max_error;
};
