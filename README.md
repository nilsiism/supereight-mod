# supereight: a fast octree library for Dense SLAM
welcome to *supereight*: a high performance template octree library and a dense
volumetric SLAM pipeline implementation.

# Related publications
This software implements the octree library and dense SLAM system presented in
our paper "Efficient Octree-Based Volumetric SLAM Supporting Signed-Distance
and Occupancy Mapping". If you publish work that relates to this software,
please cite our paper as:

@ARTICLE{VespaRAL18, 
author={E. Vespa and N. Nikolov and M. Grimm and L. Nardi and P. H. J. Kelly
and S. Leutenegger}, 
journal={IEEE Robotics and Automation Letters}, 
title={Efficient Octree-Based Volumetric SLAM Supporting Signed-Distance and
Occupancy Mapping}, year={2018}, volume={3}, number={2}, pages={1144-1151}, 
doi={10.1109/LRA.2018.2792537}, ISSN={}, month={April}}

# Licence
The core library is released under the BSD 3-clause Licence. There are part of
the this software that are released under MIT licence, see individual headers
for which licence applies.

# Project structure
supereight is made of three main different components:

* se\_core: the main header only template octree library
* se\_denseslam: the volumetric SLAM pipelines presented in [1], which can be
  compiled in a library and used in external projects.
* se\_apps: front-end applications which run the se-denseslam pipelines on
  given inputs or live camera.

Additionally, se\_tools includes the dataset generation tool and some libraries
required by se_denseslam and se_apps.

# Dependencies
The following packages are required to build the se-denseslam library:
* CMake >= 3.10
* Eigen3 
* Sophus
* TooN
* OpenMP (optional)
* GTest

The benchmarking and GUI apps additionally require:
* GLut
* OpenGL
* OpenNI2
* PAPI
* PkgConfig/Qt5
* Python/Numpy for evaluation scripts

# Build
Simply type `make` from the project root. This will create a build folder and
invoke cmake from there.
