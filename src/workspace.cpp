#include "workspace.hpp"

#include "geometry.hpp"

#include <stlplus3/filesystemSimplified/file_system.hpp>


std::string Workspace::getImagePath(int id) const
{
    return stlplus::create_filespec(root_path, view_data[id].image_path);
}

std::string Workspace::getDepthPath(int id, bool geometric) const
{
    return getPath(id, "depth", geometric);
}

std::string Workspace::getNormalPath(int id, bool geometric) const
{
    return getPath(id, "normal", geometric);
}

std::string Workspace::getPath(int id, const std::string &name, bool geometric, const std::string &extension) const
{
    static const std::string separator = "_";
    std::string postfix = geometric ? "geometric" : "photometric";
    return stlplus::create_filespec(work_path, std::to_string(id) + separator + stlplus::basename_part(view_data[id].image_path) + separator + name + separator + postfix, extension);
}

bool Workspace::getMinMaxDepth(int id, FloatT &min, FloatT &max) const
{
    const ViewData &view = view_data[id];

    min = view.min_depth;
    max = view.max_depth;

    if (min > 0 && max > min && std::isfinite(min) & std::isfinite(max))
        return true;

    std::vector<FloatT> depths;
    for (int i = 0; i < point_data.size(); ++i)
    {
        const PointData &point = point_data[i];

        if (point.projections.count(id) != 0)
        {
            const Vector3 X_= (view.R * (point.X - view.C));
            const Vector2 x = (view.K * X_).hnormalized();
            const FloatT depth = X_.z();

            if (depth > 0 && x.x() > 0 && x.x() < view.width && x.y() > 0 && x.y() < view.width)
            {
                depths.push_back(depth);
            }
        }
    }

    if (depths.empty())
        return false;

    std::sort(depths.begin(), depths.end());

    min = depths[depths.size() * 0.01];
    max = depths[depths.size() * 0.99];

    min *= 0.8;
    max *= 1.2;

    return !depths.empty();
}

bool Workspace::getSrcViewIds(int ref_id, std::vector<int> &src_view_ids) const
{
    const ViewData &ref = view_data[ref_id];

    src_view_ids = ref.src_view_ids;
    if (!src_view_ids.empty())
        return true;

    std::map<int, FloatT> score;
    for (const PointData &point : point_data)
    {
        if (point.projections.count(ref_id) == 0)
            continue;

        const Vector3 &X = point.X;
        const Vector3 &ref_C = ref.C;

        for (auto projections_iter : point.projections)
        {
            const int src_id = projections_iter.first;
            if (src_id == ref_id) continue;

            const Vector3 &src_C = view_data[src_id].C;

            const FloatT angle = angleBetweenVectors(X - ref_C, X - src_C);
            score[src_id] += std::sin(std::min(angle, static_cast<FloatT>(M_PI/2)));
        }
    }

    std::vector<std::pair<int, FloatT>> sorted(score.begin(), score.end());
    std::sort(sorted.begin(), sorted.end(), [](const std::pair<int, FloatT> &first, const std::pair<int, FloatT> &second) {
        return first.second > second.second;
    });

    for (std::pair<int, FloatT> &image_score : sorted)
        src_view_ids.push_back(image_score.first);

    return !src_view_ids.empty();
}

void Workspace::initialize()
{
    #pragma omp parallel for
    for (int view_id = 0; view_id < view_data.size(); ++view_id)
    {
        ViewData &view = view_data[view_id];

        getSrcViewIds(view_id, view.src_view_ids);
        getMinMaxDepth(view_id, view.min_depth, view.max_depth);
    }
}
