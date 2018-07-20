#include "fusion.hpp"

#include "view_utils.hpp"
#include "progress.hpp"

Fusion::Fusion(const std::shared_ptr<Workspace> &workspace)
    : workspace(workspace)
{

}

void Fusion::run(bool geometric, AbstractProgress *progress)
{
    initializeViews(geometric);

    progress->configure(views.size());

    for (int i = 0; i < views.size(); i++)
    {
        ProgressIncrementor inc(progress);

        SrcView &view = *views[i];

        if (view.depth.size() == 0 || view.normal.size() == 0 || view.color.size() == 0)
            continue;

        QueueItem item;
        item.image_id = i;

        for (int row = 0; row < view.height; row++)
        {
            for (int col = 0; col < view.width; col++)
            {
                item.x << col, row;
                fuse(item);
            }
        }
    }

    std::cout << "Exporting points" << std::endl;
    ExportPoints(points, normals, colors, stlplus::create_filespec(workspace->work_path, "fused", "ply"));
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
    if (!ref || !ref->isValid(ref_item.x)) return;

    const Vector3 X_ref = ref->unproject(ref_item.x);
    const Normal  n_ref = ref->unprojectNormal(ref_item.x);
    const auto &c_ref = ref->color(ref_item.x);

    std::vector<unsigned char> accumulated_r, accumulated_g, accumulated_b;
    std::vector<FloatT> accumulated_x,  accumulated_y,  accumulated_z;
    std::vector<FloatT>  accumulated_nx, accumulated_ny, accumulated_nz;


    std::set<std::pair<int, std::pair<int, int>>> used;
    used.insert(std::make_pair(ref_item.image_id, std::make_pair(ref_item.x(0), ref_item.x(1))));

    std::queue<QueueItem> queue;
    pupulateQueue(queue, ref_item, X_ref);

    accumulated_r.push_back(c_ref.r());
    accumulated_g.push_back(c_ref.g());
    accumulated_b.push_back(c_ref.b());
    accumulated_x.push_back(X_ref.x());
    accumulated_y.push_back(X_ref.y());
    accumulated_z.push_back(X_ref.z());
    accumulated_nx.push_back(n_ref.x());
    accumulated_ny.push_back(n_ref.y());
    accumulated_nz.push_back(n_ref.z());

    while(!queue.empty())
    {
        const Fusion::QueueItem src_item = queue.front(); queue.pop();

        std::shared_ptr<SrcView> &src = views[src_item.image_id];
        if (!src) continue;

        if (!used.insert(std::make_pair(src_item.image_id, std::make_pair(src_item.x(0), src_item.x(1)))).second)
            continue;

        if (!src->isValid(src_item.x))
            continue;

        if (!checkDepth(src->distance(X_ref), src->depth(src_item.x)))
            continue;

        const Normal n_src = src->unprojectNormal(src_item.x);

        if (!checkNormal(n_ref, n_src))
            continue;

        const Vector3 X_src = src->unproject(src_item.x);

        if (!checkReprojectionError(convertToFloat(ref_item.x), ref->project(X_src)))
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
        src->normal(src_item.x) = Normal::Zero();
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

    if (n.norm() < std::numeric_limits<FloatT>::epsilon())
        return;
    n.normalize();

    colors.push_back(c);
    points.push_back(X);
    normals.push_back(n);
}

void Fusion::pupulateQueue(std::queue<Fusion::QueueItem> &queue, const Fusion::QueueItem &item, const Vector3 &X)
{
    std::set<int> src_view_ids;
    if (!workspace->getSrcViewIds(item.image_id, src_view_ids, 20))
        return;

    for (int src_view_id : src_view_ids)
    {
        const View &view = *views[src_view_id];
        const Vector2i x = convertToInt(view.project(X));
        if (!view.isVisible(x))
            continue;

        queue.emplace(src_view_id, x, item.depth + 1);
    }
}

bool Fusion::checkDepth(const FloatT projected_depth, const FloatT src_depth) const
{
    return (std::abs(projected_depth - src_depth) / src_depth) < options.max_depth_error;
}

bool Fusion::checkNormal(const Normal &ref_normal, const Normal &src_normal) const
{
    return angleBetweenNormals(ref_normal, src_normal) < options.max_angle_error;
}

bool Fusion::checkReprojectionError(const Vector2 &ref_x, const Vector2 &ref_x_) const
{
    return (ref_x - ref_x_).norm() < options.max_reprojection_error;
}
