#pragma once

#include "workspace.hpp"

#include "stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>
#include <sstream>
#include <string>


static bool readCameras(const std::string &filename, std::map<int, Matrix3> &Ks, std::map<int, std::pair<int, int>> &sizes)
{
    std::ifstream stream(filename);
    if (!stream.is_open())
        return false;

    std::string line;
    while (std::getline(stream, line))
    {
        if (line.empty() || line.at(0) == '#')
            continue;

        std::istringstream iss(line);

        int camera_id, width, height;
        std::string camera_model;
        float fx, fy, ppx, ppy;

        iss >> camera_id >> camera_model >> width >> height >> fx >> fy >> ppx >> ppy;

        Ks[camera_id] << fx, 0, ppx, 0, fy, ppy, 0, 0, 1;
        sizes[camera_id] = std::make_pair(width, height);
    }

    return true;
}

static bool readImages(const std::string &filename, std::map<int, Matrix3> &Ks, std::map<int, std::pair<int, int>> &sizes, std::map<int, std::map<int, Vector2>> &projections, Workspace &workspace)
{
    std::ifstream stream(filename);
    if (!stream.is_open())
        return false;

    std::string line;
    while (std::getline(stream, line))
    {
        if (line.empty() || line.at(0) == '#')
            continue;

        std::istringstream iss(line);

        int image_id, camera_id;
        FloatT qw, qx, qy, qz, tx, ty, tz;
        std::string image_path;

        iss >> image_id >> qw >> qx >> qy >> qz >> tx >> ty >> tz >> camera_id >> image_path;
        Eigen::Quaternion<FloatT> q(qw, qx, qy, qz);

        image_id -= 1; // from 1-based to 0-based index

        if (workspace.view_data.size() < image_id + 1)
            workspace.view_data.resize(image_id + 1);

        ViewData &view_data = workspace.view_data[image_id];
        view_data.R = q.toRotationMatrix();
        view_data.C = -view_data.R.transpose() * Vector3(tx, ty, tz);
        view_data.K = Ks.at(camera_id);
        view_data.image_path = image_path;
        view_data.width = sizes.at(camera_id).first;
        view_data.height = sizes.at(camera_id).second;
        view_data.max_depth = std::numeric_limits<FloatT>::lowest();
        view_data.min_depth = std::numeric_limits<FloatT>::max();

        std::getline(stream, line);
        std::istringstream iss1(line);

        float x, y;
        int point_id;
        while (iss1.good())
        {
            iss1 >> x >> y >> point_id;

            if (point_id == -1) continue;

            projections[point_id][image_id] = Vector2(x, y);
        }
    }

    return true;
}

static bool readPoints(const std::string &filename, std::map<int, std::map<int, Vector2>> &projections, Workspace &workspace)
{
    std::ifstream stream(filename);
    if (!stream.is_open())
        return false;

    std::map<int, std::set<int>> src_view_ids;

    std::string line;
    while (std::getline(stream, line))
    {
        if (line.empty() || line.at(0) == '#')
            continue;

        std::istringstream iss(line);

        int point3d_id, image_id, point2d_idx;
        FloatT x, y, z, error;
        int r, g, b;

        iss >> point3d_id >> x >> y >> z >> r >> g >> b >> error;
        Vector3 X(x, y, z);


        workspace.point_data.emplace_back();
        PointData &point_data = workspace.point_data.back();
        point_data.X = X;

        std::set<int> image_ids;

        while(iss.good())
        {
            iss >> image_id >> point2d_idx;
            image_id -= 1; // from 1-based to 0-based index

            const Vector2 &projection = projections[point3d_id][image_id];
            point_data.projections[image_id] = projection;

            ViewData &view = workspace.view_data.at(image_id);
            Vector3 X_ = view.R * (X - view.C);

            view.max_depth = std::max(view.max_depth, X_.z());
            view.min_depth = std::min(view.min_depth, X_.z());

            image_ids.insert(image_id);
        }

        for (int image_id : image_ids)
        {
            //ViewData &view = workspace.view_data.at(image_id);
            src_view_ids[image_id].insert(image_ids.begin(), image_ids.end());
        }
    }

    for (int image_id = 0; image_id < workspace.view_data.size(); image_id++)
    {
        src_view_ids[image_id].erase(image_id);

        ViewData &view = workspace.view_data.at(image_id);
        view.src_view_ids.insert(view.src_view_ids.end(), src_view_ids[image_id].begin(), src_view_ids[image_id].end());
    }

    return true;
}

static bool InitializeWorkspaceEth3D(const std::string &folder, Workspace &workspace)
{
    std::map<int, Matrix3> Ks;
    std::map<int, std::pair<int, int>> sizes;
    std::map<int, std::map<int, Vector2>> projections;

    workspace.root_path = stlplus::folder_to_path(stlplus::folder_part(folder), "images");

    if (!readCameras(stlplus::create_filespec(folder, "cameras", "txt"), Ks, sizes))
        return false;
    if (!readImages(stlplus::create_filespec(folder, "images", "txt"), Ks, sizes, projections, workspace))
        return false;
    if (!readPoints(stlplus::create_filespec(folder, "points3D", "txt"), projections, workspace))
        return false;

    return true;
}

