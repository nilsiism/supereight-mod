#include <cmath>
#include "octree.hpp"
#include "math_utils.h"
#include "gtest/gtest.h"
#include "functors/axis_aligned_functor.hpp"

typedef float testT;
template <>
struct voxel_traits<testT> {
  typedef float ComputeType;
  typedef float StoredType;
  static inline ComputeType empty(){ return 0.f; }
  static inline ComputeType initValue(){ return 1.f; }
  static inline StoredType translate(const ComputeType value) {
     return value;
  }
};

float test_fun(float x, float y, float z) {
  return sq(z) + std::sin(2*x + y);
}

class InterpolationTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      unsigned size = 256;
      float dim = 5.f;
      oct_.init(size, dim); // 5 meters

      const unsigned center = 2.5f;
      const unsigned radius = center + 0.5f; 

      const float voxelsize = oct_.dim()/oct_.size();
      const float inverse_voxelsize = 1.f/voxelsize;
      unsigned leaf_level = log2(size) - log2(Octree<testT>::blockSide);
      for(float z = center - radius; z < center + radius; ++z)
        for(float y = center - radius; y < center + radius; ++y)
          for(float x = center - radius; x < center + radius; ++x) {
            const int3 vox =  make_int3(x * inverse_voxelsize, 
                y * inverse_voxelsize, z * inverse_voxelsize);
            alloc_list.push_back(oct_.hash(vox.x, vox.y, vox.z, leaf_level));
          }
      oct_.alloc_update(alloc_list.data(), alloc_list.size());
    }

  typedef Octree<testT> OctreeF;
  OctreeF oct_;
  std::vector<octlib::key_t> alloc_list;
};

TEST_F(InterpolationTest, Init) {

  auto initialise = [](auto& handler, const int3& v) {
    handler.set(test_fun(v.x, v.y, v.z));
  }; 
  iterators::functor::axis_aligned<testT, Octree, decltype(initialise)> 
    funct(oct_, initialise);
  funct.apply();

  auto test = [](auto& handler, const int3& v) {
    auto data = handler.get();
    ASSERT_EQ(data, test_fun(v.x, v.y, v.z));
  }; 
  iterators::functor::axis_aligned<testT, Octree, decltype(test)> 
    funct_test(oct_, test);
  funct_test.apply();
}

// TEST_F(InterpolationTest, InterpAtPoints) {
// 
//   auto test = [this](auto& handler, const int3& v) {
//     auto data = handler.get();
//     float interpolated = oct_.interp(make_float3(v.x, v.y, v.z), [](const auto& val){ return val; });
//     ASSERT_EQ(data, interpolated);
//   }; 
// 
//   iterators::functor::axis_aligned<testT, Octree, decltype(test)> 
//     funct_test(oct_, test);
//   funct_test.apply();
// }
