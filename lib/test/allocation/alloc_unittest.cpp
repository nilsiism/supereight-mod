#include <random>
#include "octree.hpp"
#include "math_utils.h"
#include "utils/morton_utils.hpp"
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
  const Eigen::Vector3i vox = {25, 65, 127};
  const octlib::key_t code = oct.hash(vox(0), vox(1), vox(2)); 
  octlib::key_t allocList[1] = {code};
  const float val = oct.get(vox(0), vox(1), vox(2));
  EXPECT_EQ(val, voxel_traits<float>::empty());
}

TEST(AllocationTest, SetSingleVoxel) {
  typedef Octree<float> OctreeF;
  OctreeF oct;
  oct.init(256, 5);
  const Eigen::Vector3i vox = {25, 65, 127};
  const octlib::key_t code = oct.hash(vox(0), vox(1), vox(2)); 
  octlib::key_t allocList[1] = {code};
  oct.allocate(allocList, 1);

  VoxelBlock<float> * block = oct.fetch(vox(0), vox(1), vox(2));
  float written_val = 2.f;
  block->data(vox, written_val);

  const float read_val = oct.get(vox(0), vox(1), vox(2));
  EXPECT_EQ(written_val, read_val);
}

TEST(AllocationTest, FetchOctant) {
  typedef Octree<float> OctreeF;
  OctreeF oct;
  oct.init(256, 5);
  const Eigen::Vector3i vox = {25, 65, 127};
  const uint code = oct.hash(vox(0), vox(1), vox(2)); 
  octlib::key_t allocList[1] = {code};
  oct.allocate(allocList, 1);

  const int depth = 3; /* 32 voxels per side */
  Node<float> * node = oct.fetch_octant(vox(0), vox(1), vox(2), 3);

  EXPECT_NE(node, nullptr);
}

TEST(AllocationTest, MortonPrefixMask) {

  const unsigned int max_bits = 21; 
  const unsigned int block_side = 8;
  const unsigned int size = std::pow(2, max_bits);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dis(0, size);

  constexpr int num_samples = 10;
  octlib::key_t keys[num_samples];
  octlib::key_t tempkeys[num_samples];
  Eigen::Vector3i coordinates[num_samples];

  for(int i = 0; i < num_samples; ++i) {
    const Eigen::Vector3i vox = {dis(gen), dis(gen), dis(gen)};
    coordinates[i] = Eigen::Vector3i(vox);
    const octlib::key_t code = compute_morton(vox(0), vox(1), vox(2));
    keys[i] = code;
  }

  const int max_level = log2(size);
  const int leaf_level = max_level - log2(block_side);
  const unsigned int shift = max_bits - max_level;
  int edge = size/2;
  for (int level = 0; level <= leaf_level; level++){
    const octlib::key_t mask = MASK[level + shift];
    compute_prefix(keys, tempkeys, num_samples, mask);
    for(int i = 0; i < num_samples; ++i) {
      const Eigen::Vector3i masked_vox = unpack_morton(tempkeys[i]);
      ASSERT_EQ(masked_vox(0) % edge, 0);
      ASSERT_EQ(masked_vox(1) % edge, 0);
      ASSERT_EQ(masked_vox(2) % edge, 0);
      const Eigen::Vector3i vox = coordinates[i];
      // printf("vox: %d, %d, %d\n", vox(0), vox(1), vox(2));
      // printf("masked level %d: %d, %d, %d\n", level, masked_vox(0), masked_vox(1), masked_vox(2) );
    }
    edge = edge/2;
  }
}
