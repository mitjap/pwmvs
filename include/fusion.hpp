#pragma once

#include "workspace.hpp"
#include "geometry.hpp"

#include <queue>

class AbstractProgress;

class Fusion
{
    struct QueueItem
    {
        QueueItem()
            : depth(0)
        {

        }

        QueueItem(int image_id, const Vector2i &x, int depth)
            : image_id(image_id)
            , x(x)
            , depth(depth)
        {

        }

        int image_id;
        Vector2i x;

        int depth;
    };

public:
    struct Options
    {
        FloatT max_depth_error = 0.02;
        FloatT max_angle_error = deg2rad(15);
        FloatT max_reprojection_error = 1.5;

        int min_points = 5;
    };

    Fusion(const std::shared_ptr<Workspace> &workspace);

    void run(bool geometric, AbstractProgress *progress = nullptr);

private:
    void initializeViews(bool geometric);
    void fuse(const QueueItem &ref_item);
    void pupulateQueue(std::queue<QueueItem> &queue, const QueueItem &item, const Vector3 &X);

    bool checkDepth(const FloatT projected_depth, const FloatT src_depth) const;
    bool checkNormal(const Normal &ref_normal, const Normal &src_normal) const;
    bool checkReprojectionError(const Vector2 &ref_x, const Vector2 &ref_x_) const;

public:
    std::shared_ptr<Workspace> workspace;
    Options options;

    std::vector<openMVG::image::RGBColor> colors;
    std::vector<Vector3> points;
    std::vector<Normal> normals;

    FloatT max_reprojection_error_sq;

    std::vector<unsigned char> accumulated_r, accumulated_g, accumulated_b;
    std::vector<FloatT> accumulated_x,  accumulated_y,  accumulated_z;
    std::vector<FloatT>  accumulated_nx, accumulated_ny, accumulated_nz;

private:
    std::vector<std::shared_ptr<SrcView>> views;
};
