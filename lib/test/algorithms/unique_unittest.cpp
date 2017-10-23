#include "octree.hpp"
#include "math_utils.h"
#include "utils/morton_utils.hpp"
#include "algorithms/unique.hpp"
#include "gtest/gtest.h"
#include <algorithm>

typedef float testT;
typedef unsigned int MortonType;

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

class UniqueTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      const int3 blocks[10] = {
        {56, 12, 12}, {56, 12, 15}, 
        {128, 128, 128},
        {128, 128, 125}, {128, 128, 127}, 
        {128, 136, 129}, 
        {128, 136, 127}, 
        {136, 128, 136}, 
        {128, 240, 136}, {128, 241, 136}};
      for(int i = 0; i < 10; ++i) { 
        keys[i] = oct.hash(blocks[i].x, blocks[i].y, blocks[i].z);
      }
    }
    
    MortonType keys[10];
    typedef Octree<testT> OctreeF;
    OctreeF oct;
};

class UniqueMultiscaleTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      const int3 blocks[10] = {
        {56, 12, 12}, {56, 12, 15}, 
        {128, 128, 128},
        {128, 128, 125}, {128, 128, 127}, 
        {128, 136, 129}, 
        {128, 136, 127}, 
        {136, 128, 136}, 
        {128, 240, 136}, {128, 241, 136}};
      for(int i = 0; i < 10; ++i) { 
        keys[i] = oct.hash(blocks[i].x, blocks[i].y, blocks[i].z);
      }

      keys[2] = keys[2] | 3;
      keys[3] = keys[3] | 5;
      keys[4] = keys[4] | 6;
    }
    
    MortonType keys[10];
    typedef Octree<testT> OctreeF;
    OctreeF oct;
};

TEST_F(UniqueTest, FilterDuplicates) {
  std::sort(keys, keys + 10);
  const int last = algorithms::unique(keys, 10);
  ASSERT_EQ(last, 7);
  for(int i = 1; i < last; ++i) { 
    ASSERT_TRUE(keys[i] != keys[i-1]);
  }
}

TEST_F(UniqueMultiscaleTest, FilterDuplicates) {
  std::sort(keys, keys + 10);
  /*
   * 0x1FFu extracts the last 9 bits of a morton number,
   * corresponding to the edge of a voxel block: 3*log2(VoxelBlock<T>::side)
   */
  const int last = algorithms::unique_multiscale(keys, 10, 0x1FFu);
  ASSERT_EQ(last, 7);
  for(int i = 1; i < last; ++i) { 
    std::cout << "(Key: " << (keys[i-1] & (~0x1FFu)) << ", Scale: " 
              << (keys[i-1] & 0x1FFu) << "), "; 
    ASSERT_TRUE(keys[i] != keys[i-1]);
  }
  std::cout << std::endl;
}
