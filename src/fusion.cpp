#include "fusion.hpp"

#include "view_utils.hpp"
#include "progress.hpp"

Fusion::Fusion(const std::shared_ptr<Workspace> &workspace)
    : workspace(workspace)
{
    // nothing to do
}

void Fusion::run(bool geometric, AbstractProgress *progress)
{
    max_reprojection_error_sq = options.max_reprojection_error * options.max_reprojection_error;
    cos_max_angle_error = std::cos(options.max_angle_error);

    initializeViews(geometric);

    progress->configure(views.size());


    // sort view by number of neghbours
    std::vector<int> prioritized_view_ids;
    for (int i = 0; i < views.size(); i++)
        prioritized_view_ids.push_back(i);
    std::sort(prioritized_view_ids.begin(), prioritized_view_ids.end(), [this](int a_id, int b_id) {
        return this->workspace->view_data[a_id].src_view_ids.size() > this->workspace->view_data[b_id].src_view_ids.size();
    });


    for (int src_id : prioritized_view_ids)
    {
        ProgressIncrementor inc(progress);
        std::shared_ptr<SrcView> &view = views[src_id];

        if (!view || view->depth.size() == 0 || view->normal.size() == 0 || view->color.size() == 0)
            continue;

        QueueItem item;
        item.image_id = src_id;

        for (int row = 0; row < view->height; row++)
        {
            for (int col = 0; col < view->width; col++)
            {
                item.x << col, row;
                fuse(item);
            }
        }
    }
}

void Fusion::initializeViews(bool geometric)
{
    int view_id = 0;
    for (const ViewData &view_data : workspace->view_data)
    {
        std::shared_ptr<SrcView> src = createSrcView(view_data);
        if (!loadFromFiles(*workspace, view_id, geometric, false, true, true, true, *src))
            src.reset();

        views.push_back(src);
        view_id++;
    }
}

void Fusion::fuse(const Fusion::QueueItem &ref_item)
{
    std::shared_ptr<SrcView> &ref = views[ref_item.image_id];
    if (!ref->isValid(ref_item.x)) return;

    Vector2 x = convertToFloat(ref_item.x);
    const Vector3 X_ref = ref->unproject(ref_item.x);
    const Normal  n_ref = ref->unprojectNormal(ref_item.x);
    const auto &c_ref = ref->color(ref_item.x);

    std::queue<QueueItem> queue;
    pupulateQueue(queue, ref_item, X_ref);

    accumulated_r.clear();  accumulated_r.push_back(c_ref.r());
    accumulated_g.clear();  accumulated_g.push_back(c_ref.g());
    accumulated_b.clear();  accumulated_b.push_back(c_ref.b());
    accumulated_x.clear();  accumulated_x.push_back(X_ref.x());
    accumulated_y.clear();  accumulated_y.push_back(X_ref.y());
    accumulated_z.clear();  accumulated_z.push_back(X_ref.z());
    accumulated_nx.clear(); accumulated_nx.push_back(n_ref.x());
    accumulated_ny.clear(); accumulated_ny.push_back(n_ref.y());
    accumulated_nz.clear(); accumulated_nz.push_back(n_ref.z());

    while(!queue.empty())
    {
        const Fusion::QueueItem src_item = queue.front(); queue.pop();

        std::shared_ptr<SrcView> &src = views[src_item.image_id];
        if (!src) continue;

        if (!src->isValid(src_item.x))
            continue;

        if (!checkDepth(src->distance(X_ref), src->depth(src_item.x)))
            continue;

        const Normal n_src = src->unprojectNormal(src_item.x);

        if (!checkNormal(n_ref, n_src))
            continue;

        const Vector3 X_src = src->unproject(src_item.x);

        if (!checkReprojectionError(x, ref->project(X_src)))
            continue;

        const auto &c_src = src->color(src_item.x);

        accumulated_r.push_back(c_src.r());
        accumulated_g.push_back(c_src.g());
        accumulated_b.push_back(c_src.b());
        accumulated_x.push_back(X_src.x());
        accumulated_y.push_back(X_src.y());
        accumulated_z.push_back(X_src.z());
        accumulated_nx.push_back(n_src.x());
        accumulated_ny.push_back(n_src.y());
        accumulated_nz.push_back(n_src.z());

        pupulateQueue(queue, src_item, X_src);

        // Invalidate pixel
        src->depth(src_item.x)  = 0;
        //src->normal(src_item.x) = Normal::Zero();
    }

    if (accumulated_x.size() < options.min_points)
        return;

    std::sort(accumulated_r.begin(), accumulated_r.end());
    std::sort(accumulated_g.begin(), accumulated_g.end());
    std::sort(accumulated_b.begin(), accumulated_b.end());
    std::sort(accumulated_x.begin(), accumulated_x.end());
    std::sort(accumulated_y.begin(), accumulated_y.end());
    std::sort(accumulated_z.begin(), accumulated_z.end());
    std::sort(accumulated_nx.begin(), accumulated_nx.end());
    std::sort(accumulated_ny.begin(), accumulated_ny.end());
    std::sort(accumulated_nz.begin(), accumulated_nz.end());

    int idx = accumulated_x.size() / 2;

    Vector3 X; X << accumulated_x[idx], accumulated_y[idx], accumulated_z[idx];
    Normal n; n << accumulated_nx[idx], accumulated_ny[idx], accumulated_nz[idx];
    openMVG::image::RGBColor c; c << accumulated_r[idx], accumulated_g[idx], accumulated_b[idx];

    n.normalize();

    colors.push_back(c);
    points.push_back(X);
    normals.push_back(n);
}

void Fusion::pupulateQueue(std::queue<Fusion::QueueItem> &queue, const Fusion::QueueItem &item, const Vector3 &X)
{
    for (int i = 0; i < std::min(workspace->view_data[item.image_id].src_view_ids.size(), (size_t)options.max_sources); ++i)
    {
        int src_view_id = workspace->view_data[item.image_id].src_view_ids.at(i);
        const std::shared_ptr<SrcView> &view = views[src_view_id];
        if (!view)
            continue;

        const Vector2i x = convertToInt(view->project(X));
        if (!view->isVisible(x))
            continue;

        queue.emplace(src_view_id, x, item.recursion_depth + 1);
    }
}

bool Fusion::checkDepth(const FloatT projected_depth, const FloatT src_depth) const
{
    return (std::abs(projected_depth - src_depth) / src_depth) < options.max_depth_error;
}

bool Fusion::checkNormal(const Normal &ref_normal, const Normal &src_normal) const
{
    return cosAngleBetweenNormals(ref_normal, src_normal) >= cos_max_angle_error;
}

bool Fusion::checkReprojectionError(const Vector2 &ref_x, const Vector2 &ref_x_) const
{
    return (ref_x - ref_x_).squaredNorm() < max_reprojection_error_sq;
}
