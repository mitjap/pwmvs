#include <types.hpp>
#include <image_io.hpp>

#include "stlplus3/filesystemSimplified/file_system.hpp"

#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator>

static void GetMinMaxDepth(Image<FloatT> &image, FloatT percent, FloatT &min, FloatT &max)
{
    std::vector<FloatT> data(image.data(), image.data() + (image.Width() * image.Height()));
    std::sort(data.begin(), data.end());

    int start = std::distance(data.begin(), std::find_if(data.begin(), data.end(), [&data](FloatT val) { return val != data.front(); }));
    int end = std::distance(data.rbegin(), std::find_if(data.rbegin(), data.rend(), [&data](FloatT val) { return val != data.back(); }));
    int size = data.size() - start - end;

    min = data[start + size / 2 - (size * percent)/2];
    max = data[start + size / 2 + (size * percent)/2];
}

static bool ReadColmap(const std::string &filename, Image<float> &image)
{
    std::ifstream in(filename.c_str(), std::ios::binary);
    if (!in.is_open())
    {
        std::cerr << "failed to open file: " << filename << std::endl;
        return false;
    }

    size_t width, height, depth;
    unsigned char separator;
    in >> width >> separator >> height >> separator >> depth >> separator;

    image.resize(width, height);
    in.read(reinterpret_cast<char *>(image.data()), width * height * depth * sizeof(float));

    if (in.fail())
    {
        std::cerr << "failed to read file: " << filename << std::endl;
        return false;
    }

    return true;
}

static bool ReadColmap(const std::string &filename, Image<Normal> &image)
{
    std::ifstream in(filename.c_str(), std::ios::binary);
    if (!in.is_open())
    {
        std::cerr << "failed to open file: " << filename << std::endl;
        return false;
    }

    size_t width, height, depth;
    unsigned char separator;
    in >> width >> separator >> height >> separator >> depth >> separator;

    image.resize(width, height);
    Image<FloatT> nx(width, height), ny(width, height), nz(width, height);

    in.read(reinterpret_cast<char *>(nx.data()), width * height * 1 * sizeof(float));
    in.read(reinterpret_cast<char *>(ny.data()), width * height * 1 * sizeof(float));
    in.read(reinterpret_cast<char *>(nz.data()), width * height * 1 * sizeof(float));

    if (in.fail())
    {
        std::cerr << "failed to read file: " << filename << std::endl;
        return false;
    }

    Eigen::Map<Eigen::Matrix<FloatT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic,3>>(reinterpret_cast<FloatT *>(image.data()) + 0, image.Height(), image.Width(), Eigen::Stride<Eigen::Dynamic,3>(3 * image.Width(), 3)) = nx;
    Eigen::Map<Eigen::Matrix<FloatT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic,3>>(reinterpret_cast<FloatT *>(image.data()) + 1, image.Height(), image.Width(), Eigen::Stride<Eigen::Dynamic,3>(3 * image.Width(), 3)) = ny;
    Eigen::Map<Eigen::Matrix<FloatT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic,3>>(reinterpret_cast<FloatT *>(image.data()) + 2, image.Height(), image.Width(), Eigen::Stride<Eigen::Dynamic,3>(3 * image.Width(), 3)) = nz;

    return true;
}

void InspectNormalmaps(const std::string &colmap_normal_dir)
{
    Image<Normal> normal;
    for (const std::string &file : stlplus::folder_all(colmap_normal_dir))
    {
        if (stlplus::is_folder(stlplus::folder_to_path(colmap_normal_dir, file)))
        {
            InspectNormalmaps(stlplus::folder_to_path(colmap_normal_dir, file));
            continue;
        }

        if (stlplus::extension_part(file) != std::string("bin"))
            continue;

        std::string in = stlplus::create_filespec(colmap_normal_dir, file);
        std::string out = stlplus::create_filespec(colmap_normal_dir, stlplus::basename_part(file), "jpg");
        ReadColmap(in, normal);

        DebugNormal(normal, out);
    }
}

void InspectDepthmaps(std::string colmap_depth_dir)
{
    Image<FloatT> depth;
    for (const std::string &file : stlplus::folder_all(colmap_depth_dir))
    {
        if (stlplus::is_folder(stlplus::folder_to_path(colmap_depth_dir, file)))
        {
            InspectDepthmaps(stlplus::folder_to_path(colmap_depth_dir, file));
            continue;
        }

        if (stlplus::extension_part(file) != std::string("bin"))
            continue;

        std::string in = stlplus::create_filespec(colmap_depth_dir, file);
        std::string out = stlplus::create_filespec(colmap_depth_dir, stlplus::basename_part(file), "jpg");
        ReadColmap(in, depth);

        FloatT min = 0, max = 0;
        GetMinMaxDepth(depth, 0.9, min, max);

        std::cout << "min: " << min << std::endl;
        std::cout << "max: " << max << std::endl;

        DebugImage(depth, out, min, max);
    }
}

int main()
{
    //std::string colmap_stereo_dir = "W:/20180130_Tets_270_COLMAP/Sample set #1/001_Images_GCP_file/3Dproject/temp/pc/stereo";
    //std::string colmap_stereo_dir = "W:/20180130_Tets_270_COLMAP/test_pwmvs/door_dslr_undistorted/door/colmap/stereo";
    std::string colmap_stereo_dir = "W:/20180130_Tets_270_COLMAP/test_pwmvs/piran/3Dproject/temp/pc_level_2/stereo";
    std::string colmap_depth_dir = stlplus::folder_to_path(colmap_stereo_dir, "depth_maps");
    std::string colmap_normal_dir = stlplus::folder_to_path(colmap_stereo_dir, "normal_maps");


    InspectDepthmaps(colmap_depth_dir);
    InspectNormalmaps(colmap_normal_dir);

    return EXIT_SUCCESS;

}
