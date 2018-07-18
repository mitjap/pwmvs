#pragma once

#include "types.hpp"
#include "image_io.hpp"

#include "openMVG/cameras/cameras.hpp"
#include "openMVG/cameras/Camera_undistort_image.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/image/image_resampling.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/image/image_filtering.hpp"

#include "stl/stlMap.hpp"
#include "stlplus3/filesystemSimplified/file_system.hpp"


template <typename Image>
static void Downsample(const Image& image_in,
                       const openMVG::cameras::IntrinsicBase *cam_in,
                       Image &image_out,
                       const openMVG::cameras::IntrinsicBase *cam_out)
{
    Image temp;

    openMVG::Vec kernel(5); kernel << 1, 4, 6, 4, 1; kernel /= 16;

    openMVG::image::ImageSeparableConvolution(image_in, kernel, kernel, image_out);
    openMVG::image::ImageSeparableConvolution(image_out, kernel, kernel, temp);

    TransformImage(temp, cam_in, image_out, cam_out);
}

template <typename Image>
static void Downsample(Image &in, Image &out)
{
    Image temp;
//    openMVG::image::ImageGaussianFilter(in, 0.8, temp);

    openMVG::Vec kernel(5); kernel << 1, 4, 6, 4, 1; kernel /= 16;
    std::cout << "kernel: " << kernel.transpose() << " sum: " << kernel.sum() << std::endl;
    openMVG::image::ImageSeparableConvolution(in, kernel, kernel, out);
    openMVG::image::ImageSeparableConvolution(out, kernel, kernel, temp);

    //openMVG::image::Sampler2d<openMVG::image::SamplerCubic> sampler;

//    const int new_width  = in.Width() / 2;
//    const int new_height = in.Height() / 2;

//    out.resize(new_width , new_height);

//    for (int i = 0; i < new_height; ++i)
//    {
//        for (int j = 0; j < new_width; ++j)
//        {
//            float i_ = convertToFloat(i) * 2;
//            float j_ = convertToFloat(j) * 2;

//            // offset by -0.5 because of different interpretations of floating point coordiates
//            // see https://github.com/openMVG/openMVG/issues/1271
//            out(i , j) = sampler(temp, i_ - 0.5, j_ - 0.5);
//        }
//    }

    openMVG::image::ImageDecimate(temp, out);
}

template <>
void Downsample(Image<openMVG::image::RGBColor> &image, Image<openMVG::image::RGBColor> &out)
{
    Image<unsigned char> red, green, blue;
    GetRed(image, red);
    GetGreen(image, green);
    GetBlue(image, blue);

    Downsample(red, red);
    Downsample(green, green);
    Downsample(blue, blue);

    out.resizeLike(red);

    SetRed(out, red);
    SetGreen(out, green);
    SetBlue(out, blue);
}

static std::shared_ptr<openMVG::cameras::Pinhole_Intrinsic> Downsample(const std::shared_ptr<openMVG::cameras::Pinhole_Intrinsic> &intrinsic, FloatT scale)
{
    if (auto pinhole = std::dynamic_pointer_cast<openMVG::cameras::Pinhole_Intrinsic>(intrinsic))
    {
        return std::make_shared<openMVG::cameras::Pinhole_Intrinsic>(
                    pinhole->w() * scale, pinhole->h() * scale,
                    pinhole->focal() * scale,
                    pinhole->principal_point().x() * scale, pinhole->principal_point().y() * scale);
    }
    else
    {
        std::cerr << "Unknown camera type" << std::endl;
    }

    return nullptr;
}

static std::shared_ptr<openMVG::sfm::View> Downsample(const std::shared_ptr<openMVG::sfm::View> &view, FloatT scale)
{
        std::shared_ptr<openMVG::sfm::View> downsampled = std::make_shared<openMVG::sfm::View>(*view);
        downsampled->ui_height *= scale;
        downsampled->ui_width *= scale;

        return downsampled;
}

static bool eraseObservationsWithMissingIntrinsics(openMVG::sfm::SfM_Data &sfm_data)
{
  openMVG::IndexT removed_elements = 0;

  std::set<openMVG::IndexT> intrinsic_index;
  std::transform(sfm_data.intrinsics.cbegin(), sfm_data.intrinsics.cend(),
    std::inserter(intrinsic_index, intrinsic_index.begin()), stl::RetrieveKey());

  // For each landmark:
  //  - Check if we need to keep the observations & the track
  openMVG::sfm::Landmarks::iterator itLandmarks = sfm_data.structure.begin();
  while (itLandmarks != sfm_data.structure.end())
  {
    openMVG::sfm::Observations &obs = itLandmarks->second.obs;
    openMVG::sfm::Observations::iterator itObs = obs.begin();
    while (itObs != obs.end())
    {
      const openMVG::IndexT ViewId = itObs->first;
      const std::shared_ptr<openMVG::sfm::View> &v = sfm_data.GetViews().at(ViewId);
      if (intrinsic_index.count(v->id_intrinsic) == 0)
      {
        itObs = obs.erase(itObs);
        ++removed_elements;
      }
      else
        ++itObs;
    }
    ++itLandmarks;
  }
  return removed_elements > 0;
}

static std::shared_ptr<openMVG::cameras::Pinhole_Intrinsic> ConvertToPinholeIntrinsic(const std::shared_ptr<openMVG::cameras::IntrinsicBase> &in)
{
    if (auto pinhole = std::dynamic_pointer_cast<openMVG::cameras::Pinhole_Intrinsic>(in))
    {
        return std::make_shared<openMVG::cameras::Pinhole_Intrinsic>(
                    pinhole->w(), pinhole->h(),
                    pinhole->focal(),
                    pinhole->principal_point().x(), pinhole->principal_point().y());
    }
    return nullptr;
}

static void ConvertToPinholeScene(const openMVG::sfm::SfM_Data &in, const std::string &out_root_path, openMVG::sfm::SfM_Data &out, FloatT level = 3, bool undistort = true)
{
    if (undistort)
    {
        stlplus::folder_delete(out_root_path, true);
        stlplus::folder_create(out_root_path);
    }

    out.s_root_path = out_root_path;

    // Convert intrinsics
    for (auto iter : in.intrinsics)
    {
        if (std::shared_ptr<openMVG::cameras::Pinhole_Intrinsic> pinhole = ConvertToPinholeIntrinsic(iter.second)) {
            out.intrinsics[iter.first] = Downsample(pinhole, std::pow(2, -level));

        } else {
            std::cerr << "Unable to convert intrinsic to pinhole: " << iter.first << std::endl;
        }
    }

    // Clone views
    // TODO: ATM we just copy the views since proper cloning is not supported. We could serlialize and deserialize view to get a clone.
    for (auto iter : in.views)
    {
        std::shared_ptr<openMVG::sfm::View> view = std::make_shared<openMVG::sfm::View>(*iter.second);
        out.views[iter.first] = view;
        out.views[iter.first] = Downsample(out.views[iter.first], std::pow(2, -level));
    }

    // Copy poses
    out.poses = in.poses;

    // Copy control points
    out.control_points = in.control_points;

    // Copy structure
    out.structure = in.structure;

    if (eraseObservationsWithMissingIntrinsics(out))
    {
        openMVG::sfm::eraseUnstablePosesAndObservations(out);
    }

    // //////////////// //
    // Undistort images //
    // //////////////// //
    if (undistort)
    {
        std::vector<openMVG::IndexT> keys;
        std::transform(in.views.begin(), in.views.end(), std::back_inserter(keys), [](const auto &pair){ return pair.first; });

        #pragma omp parallel for
        for (int i = 0; i < keys.size(); i++)
        {
            const std::shared_ptr<openMVG::sfm::View> &in_view = in.views.at(keys.at(i));
            const std::shared_ptr<openMVG::sfm::View> &out_view = out.views.at(keys.at(i));
            if (!in.IsPoseAndIntrinsicDefined(out_view.get()))
                continue;

            const std::string src = stlplus::create_filespec(in.s_root_path, in_view->s_Img_path);
            const std::string dst = stlplus::create_filespec(out.s_root_path, out_view->s_Img_path);

            Image<openMVG::image::RGBColor> image, temp;
            if (!Load(src, image))
            {
                std::cerr << "Unable to read image: " << src << std::endl;
                continue;
            }

            TransformImage(image, in.intrinsics.at(in_view->id_intrinsic).get(), temp, out.intrinsics.at(out_view->id_intrinsic).get(), openMVG::image::BLACK);
            image.swap(temp);

            if (!Write(dst.c_str(), image))
            {
                std::cerr << "Unable to write image: " << dst << std::endl;
                continue;
            }
        }
    }
}
