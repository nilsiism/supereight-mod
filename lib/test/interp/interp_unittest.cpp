#include "octree.hpp"
#include "math_utils.h"
#include "gtest/gtest.h"

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

class InterpTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      oct_.init(512, 5);
      const uint3 blocks[10] = {{56, 12, 254}, {87, 32, 423}, {128, 128, 128},
      {136, 128, 128}, {128, 136, 128}, {136, 136, 128}, 
      {128, 128, 136}, {136, 128, 136}, {128, 136, 136}, {136, 136, 136}};
      unsigned int alloc_list[10];
      for(int i = 0; i < 10; ++i) {
        alloc_list[i] = oct_.hash(blocks[i].x, blocks[i].y, blocks[i].z);
      }
      oct_.allocate(alloc_list, 10);
    }

  typedef Octree<testT> OctreeF;
  OctreeF oct_;
};

TEST_F(InterpTest, Init) {
  EXPECT_EQ(oct_.get(137, 138, 130), voxel_traits<testT>::initValue());
}

TEST_F(InterpTest, GatherLocal) {
  float points[8];
  VoxelBlock<testT> * block = oct_.fetch(136, 128, 136);
  const uint3 base = {136, 128, 136};
  gather_local(block, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(InterpTest, Gather4Case1) {
  float points[8];
  const uint3 base = {135, 128, 136};
  VoxelBlock<testT> * block = oct_.fetch(base.x, base.y, base.z);

  const unsigned int offs1[4] = {0, 2, 4, 6};
  gather_4(block, base, [](const auto& val){ return val; }, offs1, points);

  const unsigned int offs2[4] = {1, 3, 5, 7};
  const uint3 base1 = base + offs2[0];
  block = oct_.fetch(base1.x, base1.y, base1.z);
  gather_4(block, base, [](const auto& val){ return val; }, offs2, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}
