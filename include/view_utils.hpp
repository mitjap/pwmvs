#pragma once

#include "types.hpp"
#include "workspace.hpp"
#include "image_io.hpp"

#include <stlplus3/filesystemSimplified/file_system.hpp>


static RefView *createRefView(const ViewData &ref_view)
{
    RefView *ref = new RefView(ref_view.width, ref_view.height);
    ref->K = ref_view.K;
    ref->K_inv = ref_view.K.inverse();

    return ref;
}

static SrcView *createSrcView(const ViewData &src_view)
{
    SrcView *src = new SrcView(src_view.width, src_view.height);

    src->K = src_view.K;
    src->K_inv = src_view.K.inverse();

    // Relative pose
    src->R = src_view.R;
    src->t = src_view.t();
    src->initialize();

    return src;
}

static SrcView *createSrcView(const ViewData &ref_view, const ViewData &src_view)
{
    SrcView *src = createSrcView(src_view);

    // Relative pose
    src->R = src_view.R * ref_view.R.transpose();
    src->C = ref_view.t() + ref_view.R * src_view.C;
    src->t = -src->R * src->C;
    src->initialize();

    return src;
}

static bool isComplete(const Workspace &workspace, int view_id, bool geometric)
{
    return stlplus::file_exists(workspace.getDepthPath(view_id, geometric))
            && stlplus::file_exists(workspace.getNormalPath(view_id, geometric));
}

static bool loadFromFiles(const Workspace &workspace, int view_id, bool geometric, bool image, bool color, bool depths, bool normals, View &view)
{
    if (image && !Read(workspace.getImagePath(view_id), view.image))
        return false;
    if (color && !Read(workspace.getImagePath(view_id), view.color))
        return false;
    if (depths && !Read(workspace.getDepthPath(view_id, geometric), view.depth))
        return false;
    if (normals && !Read(workspace.getNormalPath(view_id, geometric), view.normal))
        return false;

    return true;
}

static bool saveToFiles(const Workspace &workspace, int view_id, bool depths, bool normals, bool geometric, const View &view)
{
    if (depths && !Write(workspace.getDepthPath(view_id, geometric), view.depth))
        return false;
    if (normals && !Write(workspace.getNormalPath(view_id, geometric), view.normal))
        return false;

    return true;
}

static bool debugToFiles(const Workspace &workspace, int view_id, bool depths, bool normals, bool geometric, const View &view, FloatT min_depth = 0, FloatT max_depth = 0)
{
    if (min_depth == 0) min_depth = workspace.view_data[view_id].min_depth;
    if (max_depth == 0) max_depth = workspace.view_data[view_id].max_depth;

    const ViewData view_ = workspace.view_data[view_id];

    ExportPoints(view, workspace.getPath(view_id, "cloud", geometric, "ply"), view_.R.transpose(), -view_.R.transpose() * view_.t());

    if (depths)
        DebugImage(view.depth, workspace.getPath(view_id, "depth_debug", geometric, "jpg"), min_depth, max_depth);
    if (normals)
        DebugNormal(view.normal, workspace.getPath(view_id, "normal_debug", geometric, "jpg"));

    return true;
}
