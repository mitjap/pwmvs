cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(pwmvs VERSION 0.0.1 LANGUAGES CXX)

find_package(OpenMVG REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)

add_library(pwmvs
    src/fusion.cpp
    src/pwmvs.cpp
    src/pwmvs_controller.cpp
    src/workspace.cpp
    src/progress.cpp
)

target_include_directories(pwmvs
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        src
)

target_link_libraries(pwmvs
    PUBLIC
        OpenMVG::openMVG_multiview
        OpenMVG::openMVG_sfm
        OpenMVG::openMVG_image
        OpenMVG::openMVG_features
        Eigen3::Eigen
    PRIVATE
        OpenMP::OpenMP_CXX
)
target_compile_features(pwmvs
    PUBLIC
        cxx_generic_lambdas
)


add_executable(run "src/main.cpp")
target_include_directories(run PRIVATE src)
target_link_libraries(run PRIVATE pwmvs)

add_executable(colmap_investigator "src/colmap_investigator.cpp")
target_include_directories(colmap_investigator PRIVATE src)
target_link_libraries(colmap_investigator PRIVATE pwmvs)

install(DIRECTORY include/ DESTINATION include)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/core_export.h DESTINATION include)


install(
    TARGETS pwmvs
    EXPORT "pwmvsTargets"
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)
install(EXPORT "pwmvsTargets"
    DESTINATION share/pwmvs/cmake
)

# install config files
include(CMakePackageConfigHelpers)
configure_package_config_file(pwmvsConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/pwmvsConfig.cmake
  INSTALL_DESTINATION share/pwmvs/cmake
)
write_basic_package_version_file(
  pwmvsConfigVersion.cmake
  VERSION 0.0.1
  COMPATIBILITY SameMajorVersion
)
install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/pwmvsConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/pwmvsConfigVersion.cmake
    DESTINATION share/pwmvs/cmake
)
