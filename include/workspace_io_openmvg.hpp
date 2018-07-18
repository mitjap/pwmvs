#pragma once

#include "types.hpp"
#include "workspace.hpp"
#include "sfm_data_undistorter.hpp"


#include <openMVG/image/image_io.hpp>
#include <openMVG/cameras/cameras.hpp>
#include <openMVG/cameras/Camera_Pinhole.hpp>
#include <openMVG/sfm/sfm.hpp>

static void InitializeWorkspace(openMVG::sfm::SfM_Data &sfm_data, Workspace &workspace)
{
    workspace.root_path = sfm_data.s_root_path;

    std::map<openMVG::IndexT, int> id_map;
    for (auto openmvg_view_iter : sfm_data.views)
    {
        const std::shared_ptr<openMVG::sfm::View> &openmvg_view = openmvg_view_iter.second;
        if (!sfm_data.IsPoseAndIntrinsicDefined(openmvg_view.get()))
            continue;

        ViewData view;
        view.image_path = openmvg_view->s_Img_path;
        view.height = openmvg_view->ui_height;
        view.width = openmvg_view->ui_width;

        openMVG::geometry::Pose3 openmvg_pose = sfm_data.poses[openmvg_view->id_pose];
        view.R = openmvg_pose.rotation().cast<FloatT>();
        view.C = openmvg_pose.center().cast<FloatT>();

        if (std::shared_ptr<openMVG::cameras::Pinhole_Intrinsic> openmvg_pinhole = std::dynamic_pointer_cast<openMVG::cameras::Pinhole_Intrinsic>(sfm_data.intrinsics[openmvg_view->id_intrinsic]))
        {
            if (openmvg_pinhole->have_disto())
                continue;

            view.K = openmvg_pinhole->K().cast<FloatT>();
        }
        else
        {
            continue;
        }

        id_map[openmvg_view->id_view] = workspace.view_data.size();
        workspace.view_data.push_back(view);
    }

    for (auto openmvg_landmark_iter : sfm_data.structure)
    {
        const openMVG::sfm::Landmark &landmark = openmvg_landmark_iter.second;

        PointData point_data;
        point_data.X = landmark.X.cast<FloatT>();

        for (auto openmvg_observation_iter : landmark.obs)
        {
            const openMVG::IndexT id_view = openmvg_observation_iter.first;
            const openMVG::sfm::Observation &openmvg_observation = openmvg_observation_iter.second;

            if (id_map.count(id_view) > 0)
            {
                point_data.projections[id_map.at(id_view)] = Vector2(convertFromOpenMVG(openmvg_observation.x.x()), convertFromOpenMVG(openmvg_observation.x.y()));
            }
        }

        if (!point_data.projections.empty())
            workspace.point_data.push_back(point_data);
    }
}

static void InitializeWorkspaceOpenMVG(const std::string &sfm_data_path, const std::string &root_path, Workspace &workspace, FloatT level)
{
    openMVG::sfm::SfM_Data sfm_data_in, sfm_data;
    openMVG::sfm::Load(sfm_data_in, sfm_data_path, openMVG::sfm::ESfM_Data::ALL);

    if (stlplus::is_relative_path(sfm_data_in.s_root_path))
        sfm_data_in.s_root_path = stlplus::filespec_to_path(stlplus::folder_part(sfm_data_path), sfm_data_in.s_root_path);

    stlplus::folder_create(root_path);

    ConvertToPinholeScene(sfm_data_in, root_path, sfm_data, level, true);
    InitializeWorkspace(sfm_data, workspace);
}
