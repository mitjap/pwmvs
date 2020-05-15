# PWMVS
CPU implementation of "Pixelwise View Selection for Unstructured Multi-View Stereo (Schönberger et al.)" 

For GPU implementation go to https://github.com/colmap/colmap.

# How to build
## Dependencies
- [OpenMVG](https://github.com/openMVG/openMVG)
- [Eigen3](eigen.tuxfamily.org)
- [OpenMP](https://www.openmp.org/)

First make sure you have all dependencies installed on your system. 

`sudo apt-get install libomp-dev`

`sudo apt-get install libeigen3-dev`

To build OpenMVG follow [these](https://github.com/openMVG/openMVG/blob/develop/BUILD.md) instructions. I will refer to OpenMVG installation directory as `<OpenMVG_install_dir>`.

```
git clone git@github.com:mitjap/pwmvs.git

mkdir pwmvs-build && cd pwmvs-build

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../ pwmvs-install -DOpenMVG_DIR=<OpenMVG_install_dir>/share/openMVG/cmake ../pwmvs

cmake --build . --target install
```


# How to include PWMVS to your software
Including PWMVS in your software should be fairly simple. I'll refer to PWMVS installation directory as `<PWMVS_install_dir>`

Just add this to your `CMakeLists.txt` file:
```
find_package(pwmvs)
target_link_libraries(<your_target> PRIVATE pwmvs)
```
You will also need to specify PWMVS installation directory with:
```
pwmvs_DIR=<PWMVS_install_dir>/share/pwmvs/cmake
```

# How to use PWMVS
For complete example take a look at `main.cpp`. To run dense reconstruction on your OpenMVG project export `sfm_data.json` in same folder as images, set proper path in `main.cpp` and run. 


