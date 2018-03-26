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
  static inline ComputeType initValue(){ return 0.f; }
  static inline StoredType translate(const ComputeType value) {
     return value;
  }
};

class AxisAlignedTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      unsigned size = 256;
      float dim = 5.f;
      oct_.init(size, dim); // 5 meters

      const unsigned center = 2.5f;

      const float voxelsize = oct_.dim()/oct_.size();
      const float inverse_voxelsize = 1.f/voxelsize;
      const int band = 1 * inverse_voxelsize;
      const int3 offset = make_int3(oct_.size()/2 - band/2);
      unsigned leaf_level = log2(size) - log2(Octree<testT>::blockSide);
      for(int z = 0; z < band; ++z) {
        for(int y = 0; y < band; ++y) {
          for(int x = 0; x < band; ++x) {
            const int3 vox =  make_int3(x + offset.x, y + offset.y, z + offset.z);
            alloc_list.push_back(oct_.hash(vox.x, vox.y, vox.z, leaf_level));
          }
        }
      }
      oct_.alloc_update(alloc_list.data(), alloc_list.size());
    }

  typedef Octree<testT> OctreeF;
  OctreeF oct_;
  std::vector<octlib::key_t> alloc_list;
};

TEST_F(AxisAlignedTest, Init) {

  auto initialise = [](auto& handler, const int3&) {
    handler.set(1.f);
  }; 

  auto test = [](auto& handler, const int3&) {
    auto data = handler.get();
    ASSERT_EQ(data, 1.f);
  }; 

  iterators::functor::axis_aligned<testT, Octree, decltype(initialise)> 
    funct(oct_, initialise);
  funct.apply();

  iterators::functor::axis_aligned<testT, Octree, decltype(test)> 
    funct_test(oct_, test);
  funct_test.apply();
}

TEST_F(AxisAlignedTest, BBoxTest) {

  auto set_to_ten = [](auto& handler, const int3& v) {
        if(v.x > 100 && v.x < 150 && v.y > 100 && v.y < 150 && v.z > 100 
            && v.z < 150) {
          handler.set(10.f);
        }
    };

  iterators::functor::axis_aligned<testT, Octree, decltype(set_to_ten)> 
    funct(oct_, set_to_ten);
  funct.apply();

  for(int z = 50; z < 200; ++z)
    for(int y = 50; y < 200; ++y)
      for(int x = 50; x < 200; ++x) {
        auto data = oct_.get(x, y, z);
        if(x > 100 && x < 150 && y > 100 && y < 150 && z > 100 
            && z < 150) {
          if(data != 10.f) {
            oct_.get(x, y, z);
          }
        }
      }

}
