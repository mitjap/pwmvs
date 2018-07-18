#pragma once

#include "types.hpp"

#include <string>
#include <set>

struct ViewData
{
    std::string image_path;

    int width, height;

    Matrix3 K;
    Matrix3 R;
    Vector3 C;

    FloatT min_depth = 0;
    FloatT max_depth = 0;

    std::set<int> src_view_ids;

    Vector3 t() const { return -R * C; }
};

struct PointData
{
    Vector3 X;
    std::map<int, Vector2> projections;
};

struct Workspace
{
    std::string getImagePath(int id) const;
    std::string getDepthPath(int id, bool geometric) const;
    std::string getNormalPath(int id, bool geometric) const;
    std::string getPath(int id, const std::string &name, bool geometric, const std::string &extension = "bin") const;

    std::vector<ViewData> view_data;
    std::vector<PointData> point_data;

    bool getMinMaxDepth(int id, FloatT &min, FloatT &max) const;
    bool getSrcViewIds(int ref_id, std::set<int> &src_view_ids, int max_source_views) const;

    void initialize(int max_source_views);

    std::string root_path;
    std::string work_path;
};
