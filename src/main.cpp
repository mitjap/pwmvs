#include "pwmvs_controller.hpp"
#include "pwmvs.hpp"
#include "fusion.hpp"
#include "workspace_io_openmvg.hpp"
#include "workspace_io_eth3d.hpp"

#include <openMVG/sfm/sfm.hpp>

#include "progress.hpp"

int main()
{
    //const std::string eth3d_path = "W:/20180130_Tets_270_COLMAP/test_pwmvs/door_dslr_undistorted/door/dslr_calibration_undistorted";
    //const std::string eth3d_path = "W:/20180130_Tets_270_COLMAP/test_pwmvs/delivery_area_rig_undistorted/delivery_area/rig_calibration_undistorted";
    //InitializeWorkspaceEth3D(eth3d_path, pwmvs_controller.workspace);

    //const std::string sfm_data_path = "W:/20180130_Tets_270_COLMAP/test_pwmvs/policija/cameras.json";
    const std::string sfm_data_path = "W:/20180130_Tets_270_COLMAP/test_pwmvs/moseja/cameras.json";
    //const std::string sfm_data_path = "M:/DEVELOPMENT/pwmvs_test/piran/cameras.json";

    std::shared_ptr<Workspace> workspace = std::make_shared<Workspace>();

    const std::string root_path = stlplus::folder_to_path(stlplus::folder_part(sfm_data_path), "undistorted");
    InitializeWorkspaceOpenMVG(sfm_data_path, root_path, *workspace, 0.5);

    workspace->work_path = workspace->root_path + "/pwmvs";
    workspace->initialize(20);

    bool geometric = false;
    ConsoleProgress pwmvs_progress;

    Controller<PWMVS> pwmvs_controller(workspace);
    pwmvs_controller.pwmvs_options.monte_carlo_samples = 15;
    pwmvs_controller.pwmvs_options.num_iterations = 5;
    pwmvs_controller.pwmvs_options.window_size = 5;
    pwmvs_controller.pwmvs_options.perturbation = 1;
    pwmvs_controller.pwmvs_options.filter_photometric_consistency = true;
    pwmvs_controller.pwmvs_options.filter_min_color_similarity = 0.2;
    pwmvs_controller.pwmvs_options.filter_min_num_consistent = 2;
    pwmvs_controller.pwmvs_options.filter_geometric_consistency = true;

    if (!pwmvs_controller.run(geometric, &pwmvs_progress))
        return EXIT_FAILURE;

    Fusion fusion(workspace);
    fusion.options.min_points = 4;
    fusion.options.max_reprojection_error = 1.3;
    fusion.run(geometric, &pwmvs_progress);


    return EXIT_SUCCESS;
}
