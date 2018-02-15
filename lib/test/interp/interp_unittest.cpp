#include "octree.hpp"
#include "math_utils.h"
#include "interp_gather.hpp"
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
      const int3 blocks[10] = {{56, 12, 254}, {87, 32, 423}, {128, 128, 128},
      {136, 128, 128}, {128, 136, 128}, {136, 136, 128}, 
      {128, 128, 136}, {136, 128, 136}, {128, 136, 136}, {136, 136, 136}};
      morton_type alloc_list[10];
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
  const int3 base = {136, 128, 136};
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(InterpTest, ZCrosses) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const int3 base = {132, 128, 135};
  unsigned int crossmask = ((base.x % blockSize) == blockSize - 1 << 2) | 
                           ((base.y % blockSize) == blockSize - 1 << 1) |
                            (base.z % blockSize) == blockSize - 1;
  ASSERT_EQ(crossmask, 1);
  VoxelBlock<testT> * block = oct_.fetch(base.x, base.y, base.z);
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(InterpTest, YCrosses) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const int3 base = {132, 135, 132};
  unsigned int crossmask = ((base.x % blockSize == blockSize - 1) << 2) | 
                           ((base.y % blockSize == blockSize - 1) << 1) |
                            ((base.z % blockSize) == blockSize - 1);
  ASSERT_EQ(crossmask, 2);
  VoxelBlock<testT> * block = oct_.fetch(base.x, base.y, base.z);
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(InterpTest, XCrosses) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const int3 base = {135, 132, 132};
  unsigned int crossmask = ((base.x % blockSize == blockSize - 1) << 2) | 
                           ((base.y % blockSize == blockSize - 1) << 1) |
                            ((base.z % blockSize) == blockSize - 1);
  ASSERT_EQ(crossmask, 4);
  VoxelBlock<testT> * block = oct_.fetch(base.x, base.y, base.z);
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(InterpTest, YZCross) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const int3 base = {129, 135, 135};
  unsigned int crossmask = ((base.x % blockSize == blockSize - 1) << 2) | 
                           ((base.y % blockSize == blockSize - 1) << 1) |
                            ((base.z % blockSize) == blockSize - 1);
  ASSERT_EQ(crossmask, 3);
  VoxelBlock<testT> * block = oct_.fetch(base.x, base.y, base.z);
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(InterpTest, XZCross) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const int3 base = {135, 131, 135};
  unsigned int crossmask = ((base.x % blockSize == blockSize - 1) << 2) | 
                           ((base.y % blockSize == blockSize - 1) << 1) |
                            ((base.z % blockSize) == blockSize - 1);
  ASSERT_EQ(crossmask, 5);
  VoxelBlock<testT> * block = oct_.fetch(base.x, base.y, base.z);
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(InterpTest, XYCross) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const int3 base = {135, 135, 138};
  unsigned int crossmask = ((base.x % blockSize == blockSize - 1) << 2) | 
                           ((base.y % blockSize == blockSize - 1) << 1) |
                            ((base.z % blockSize) == blockSize - 1);
  ASSERT_EQ(crossmask, 6);
  VoxelBlock<testT> * block = oct_.fetch(base.x, base.y, base.z);
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(InterpTest, AllCross) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const int3 base = {135, 135, 135};
  unsigned int crossmask = ((base.x % blockSize == blockSize - 1) << 2) | 
                           ((base.y % blockSize == blockSize - 1) << 1) |
                            ((base.z % blockSize) == blockSize - 1);
  ASSERT_EQ(crossmask, 7);
  VoxelBlock<testT> * block = oct_.fetch(base.x, base.y, base.z);
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}
