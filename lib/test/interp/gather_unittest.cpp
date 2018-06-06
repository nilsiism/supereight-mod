#include "octree.hpp"
#include "utils/se_common.h"
#include "interpolation/interp_gather.hpp"
#include "gtest/gtest.h"

typedef float testT;

template <>
struct voxel_traits<testT> {
  typedef float value_type;
  static inline value_type empty(){ return 0.f; }
  static inline value_type initValue(){ return 1.f; }
};

class GatherTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      oct_.init(512, 5);
      const Eigen::Vector3i blocks[10] = {{56, 12, 254}, {87, 32, 423}, {128, 128, 128},
      {136, 128, 128}, {128, 136, 128}, {136, 136, 128}, 
      {128, 128, 136}, {136, 128, 136}, {128, 136, 136}, {136, 136, 136}};
      se::key_t alloc_list[10];
      for(int i = 0; i < 10; ++i) {
        alloc_list[i] = oct_.hash(blocks[i](0), blocks[i](1), blocks[i](2));
      }
      oct_.allocate(alloc_list, 10);
    }

  typedef Octree<testT> OctreeF;
  OctreeF oct_;
};

TEST_F(GatherTest, Init) {
  EXPECT_EQ(oct_.get(137, 138, 130), voxel_traits<testT>::initValue());
}

TEST_F(GatherTest, GatherLocal) {
  float points[8];
  const Eigen::Vector3i base = {136, 128, 136};
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(GatherTest, ZCrosses) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const Eigen::Vector3i base = {132, 128, 135};
  unsigned int crossmask = ((base(0) % blockSize) == blockSize - 1 << 2) | 
                           ((base(1) % blockSize) == blockSize - 1 << 1) |
                            (base(2) % blockSize) == blockSize - 1;
  ASSERT_EQ(crossmask, 1);
  VoxelBlock<testT> * block = oct_.fetch(base(0), base(1), base(2));
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(GatherTest, YCrosses) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const Eigen::Vector3i base = {132, 135, 132};
  unsigned int crossmask = ((base(0) % blockSize == blockSize - 1) << 2) | 
                           ((base(1) % blockSize == blockSize - 1) << 1) |
                            ((base(2) % blockSize) == blockSize - 1);
  ASSERT_EQ(crossmask, 2);
  VoxelBlock<testT> * block = oct_.fetch(base(0), base(1), base(2));
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(GatherTest, XCrosses) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const Eigen::Vector3i base = {135, 132, 132};
  unsigned int crossmask = ((base(0) % blockSize == blockSize - 1) << 2) | 
                           ((base(1) % blockSize == blockSize - 1) << 1) |
                            ((base(2) % blockSize) == blockSize - 1);
  ASSERT_EQ(crossmask, 4);
  VoxelBlock<testT> * block = oct_.fetch(base(0), base(1), base(2));
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(GatherTest, YZCross) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const Eigen::Vector3i base = {129, 135, 135};
  unsigned int crossmask = ((base(0) % blockSize == blockSize - 1) << 2) | 
                           ((base(1) % blockSize == blockSize - 1) << 1) |
                            ((base(2) % blockSize) == blockSize - 1);
  ASSERT_EQ(crossmask, 3);
  VoxelBlock<testT> * block = oct_.fetch(base(0), base(1), base(2));
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(GatherTest, XZCross) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const Eigen::Vector3i base = {135, 131, 135};
  unsigned int crossmask = ((base(0) % blockSize == blockSize - 1) << 2) | 
                           ((base(1) % blockSize == blockSize - 1) << 1) |
                            ((base(2) % blockSize) == blockSize - 1);
  ASSERT_EQ(crossmask, 5);
  VoxelBlock<testT> * block = oct_.fetch(base(0), base(1), base(2));
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(GatherTest, XYCross) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const Eigen::Vector3i base = {135, 135, 138};
  unsigned int crossmask = ((base(0) % blockSize == blockSize - 1) << 2) | 
                           ((base(1) % blockSize == blockSize - 1) << 1) |
                            ((base(2) % blockSize) == blockSize - 1);
  ASSERT_EQ(crossmask, 6);
  VoxelBlock<testT> * block = oct_.fetch(base(0), base(1), base(2));
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}

TEST_F(GatherTest, AllCross) {
  float points[8];
  const uint blockSize = VoxelBlock<testT>::side;
  const Eigen::Vector3i base = {135, 135, 135};
  unsigned int crossmask = ((base(0) % blockSize == blockSize - 1) << 2) | 
                           ((base(1) % blockSize == blockSize - 1) << 1) |
                            ((base(2) % blockSize) == blockSize - 1);
  ASSERT_EQ(crossmask, 7);
  VoxelBlock<testT> * block = oct_.fetch(base(0), base(1), base(2));
  gather_points(oct_, base, [](const auto& val){ return val; }, points);

  for(int i = 0; i < 8; ++i) {
    EXPECT_EQ(points[i], voxel_traits<testT>::initValue());
  }
}
