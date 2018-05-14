#include <cmath>
#include "octree.hpp"
#include "utils/eigen_helper.h"
#include "gtest/gtest.h"
#include "functors/axis_aligned_functor.hpp"

typedef float testT;
template <>
struct voxel_traits<testT> {
  typedef float1 ComputeType;
  typedef float1 StoredType;
  static inline ComputeType empty(){ return {0.f}; }
  static inline ComputeType initValue(){ return {1.f}; }
  static inline StoredType translate(const ComputeType value) {
     return value;
  }
};

float test_fun(float x, float y, float z) {
  return octlib::math::sq(z) + std::sin(2*x + y);
}

class InterpolationTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      unsigned size = 512;
      float dim = 5.f;
      oct_.init(size, dim); // 5 meters

      const unsigned center = 2.5f;
      const unsigned radius = center + 0.5f; 

      const float voxelsize = oct_.dim()/oct_.size();
      const float inverse_voxelsize = 1.f/voxelsize;
      const int band = 1 * inverse_voxelsize;
      const Eigen::Vector3i offset = 
        Eigen::Vector3i::Constant(oct_.size()/2 - band/2);
      unsigned leaf_level = log2(size) - log2(Octree<testT>::blockSide);
      for(int z = 0; z < band; ++z) {
        for(int y = 0; y < band; ++y) {
          for(int x = 0; x < band; ++x) {
            const Eigen::Vector3i vox =  Eigen::Vector3i(x + offset(0), 
                y + offset(1), z + offset(2));
            alloc_list.push_back(oct_.hash(vox(0), vox(1), vox(2), leaf_level));
          }
        }
      }
      oct_.alloc_update(alloc_list.data(), alloc_list.size());
    }

  typedef Octree<testT> OctreeF;
  OctreeF oct_;
  std::vector<octlib::key_t> alloc_list;
};

TEST_F(InterpolationTest, Init) {

  auto initialise = [](auto& handler, const Eigen::Vector3i& v) {
    float1 data;
    data.x = test_fun(v(0), v(1), v(2));
    handler.set(data);
  }; 

  iterators::functor::axis_aligned<testT, Octree, decltype(initialise)> 
    funct(oct_, initialise);
  funct.apply();

  auto test = [](auto& handler, const Eigen::Vector3i& v) {
    auto data = handler.get();
    ASSERT_EQ(data.x, test_fun(v(0), v(1), v(2)));
  }; 
  iterators::functor::axis_aligned<testT, Octree, decltype(test)> 
    funct_test(oct_, test);
  funct_test.apply();

  // std::stringstream f;
  // f << "./analytical_function.vtk";
  // save3DSlice(oct_, Eigen::Vector3i(0, oct_.size()/2, 0),
  //     Eigen::Vector3i(oct_.size(), oct_.size()/2 + 1, oct_.size()), 
  //     Eigen::Vector3i(oct_.size()), f.str().c_str());
  // f.str("");
  // f.clear();
}

// TEST_F(InterpolationTest, InterpAtPoints) {
// 
//   auto test = [this](auto& handler, const Eigen::Vector3i& v) {
//     auto data = handler.get();
//     float interpolated = oct_.interp(make_float3(v(0), v(1), v(2)), [](const auto& val){ return val(0); });
//     ASSERT_EQ(data(0), interpolated);
//   }; 
// 
//   iterators::functor::axis_aligned<testT, Octree, decltype(test)> 
//     funct_test(oct_, test);
//   funct_test.apply();
// }
