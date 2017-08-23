#include "octree.hpp"
#include "math_utils.h"
#include "gtest/gtest.h"

template <>
struct voxel_traits<float> {
  typedef float ComputeType;
  typedef float StoredType;
  static inline ComputeType empty(){ return 0.f; }
  static inline ComputeType initValue(){ return 0.f; }
  static inline StoredType translate(const ComputeType value) {
     return value;
  }
};

TEST(AllocationTest, EmptySingleVoxel) {
  typedef Octree<float> OctreeF;
  OctreeF oct;
  oct.init(256, 5);
  const int3 vox = {25, 65, 127};
  const morton_type code = oct.hash(vox.x, vox.y, vox.z); 
  morton_type allocList[1] = {code};
  const float val = oct.get(vox.x, vox.y, vox.z);
  EXPECT_EQ(val, voxel_traits<float>::empty());
}

TEST(AllocationTest, SetSingleVoxel) {
  typedef Octree<float> OctreeF;
  OctreeF oct;
  oct.init(256, 5);
  const int3 vox = {25, 65, 127};
  const morton_type code = oct.hash(vox.x, vox.y, vox.z); 
  morton_type allocList[1] = {code};
  oct.allocate(allocList, 1);

  VoxelBlock<float> * block = oct.fetch(vox.x, vox.y, vox.z);
  float written_val = 2.f;
  block->data(vox, written_val);

  const float read_val = oct.get(vox.x, vox.y, vox.z);
  EXPECT_EQ(written_val, read_val);
}

TEST(AllocationTest4096, SetSingleVoxel) {
  typedef Octree<float> OctreeF;
  OctreeF oct;
  oct.init(4096, 5);
  const int3 vox = {25, 2022, 1980};
  const morton_type code = oct.hash(vox.x, vox.y, vox.z); 
  morton_type allocList[1] = {code};
  oct.allocate(allocList, 1);

  VoxelBlock<float> * block = oct.fetch(vox.x, vox.y, vox.z);
  ASSERT_TRUE(block);
  float written_val = 2.f;
  block->data(vox, written_val);

  const float read_val = oct.get(vox.x, vox.y, vox.z);
  EXPECT_EQ(written_val, read_val);
}
