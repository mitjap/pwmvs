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
        FloatT max_depth_error = 0.01;
        FloatT max_angle_error = deg2rad(10);
        FloatT max_reprojection_error = 2;

        int min_points = 5;
    };

    Fusion(const std::shared_ptr<Workspace> &workspace);
    virtual ~Fusion();

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

private:
    std::vector<SrcView *> views;
};
