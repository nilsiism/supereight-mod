#include <fstream>
#include <string>
#include <random>
#include "se_serialise.hpp"
#include "node.hpp"
#include "octree.hpp"
#include "gtest/gtest.h"

typedef float testT;
template <>
struct voxel_traits<testT> {
  typedef float value_type;
  static inline value_type empty(){ return 0.f; }
  static inline value_type initValue(){ return 10.f; }
};

class Occupancy;
template <>
struct voxel_traits<Occupancy> {
  typedef struct { 
    float x;
    double y;
  } value_type;
  static inline value_type empty(){ return {0.f, 0.}; }
  static inline value_type initValue(){ return {1.f, 0.}; }
};

TEST(SerialiseUnitTest, WriteReadNode) {
  std::string filename = "test.bin";
  {
    std::ofstream os (filename, std::ios::binary); 
    se::Node<testT> octant;
    octant.code = 24;
    octant.side = 256;
    for(int i = 0; i < 8; ++i)
      octant.value_[i] =  5.f;
    se::internal::serialise(os, octant);
  }

  {
    std::ifstream is(filename, std::ios::binary);
    se::Node<testT> octant;
    se::internal::deserialise(octant, is);
    ASSERT_EQ(octant.code, 24);
    ASSERT_EQ(octant.side, 256);
    for(int i = 0; i < 8; ++i)
      ASSERT_EQ(octant.value_[i], 5.f);
  }
}

TEST(SerialiseUnitTest, WriteReadBlock) {
  std::string filename = "test.bin";
  {
    std::ofstream os (filename, std::ios::binary); 
    se::VoxelBlock<testT> octant;
    octant.code = 24;
    octant.coordinates(Eigen::Vector3i(40, 48, 52));
    for(int i = 0; i < 512; ++i)
      octant.data(i, 5.f);
    se::internal::serialise(os, octant);
  }

  {
    std::ifstream is(filename, std::ios::binary);
    se::VoxelBlock<testT> octant;
    se::internal::deserialise(octant, is);
    ASSERT_EQ(octant.code, 24);
    ASSERT_TRUE(octant.coordinates() == Eigen::Vector3i(40, 48, 52));
    for(int i = 0; i < 512; ++i)
      ASSERT_EQ(octant.data(i), 5.f);
  }
}

TEST(SerialiseUnitTest, WriteReadBlockStruct) {
  std::string filename = "test.bin";
  {
    std::ofstream os (filename, std::ios::binary); 
    se::VoxelBlock<Occupancy> octant;
    octant.code = 24;
    octant.coordinates(Eigen::Vector3i(40, 48, 52));
    for(int i = 0; i < 512; ++i)
      octant.data(i, {5.f, 2.});
    se::internal::serialise(os, octant);
  }

  {
    std::ifstream is(filename, std::ios::binary);
    se::VoxelBlock<Occupancy> octant;
    se::internal::deserialise(octant, is);
    ASSERT_EQ(octant.code, 24);
    ASSERT_TRUE(octant.coordinates() == Eigen::Vector3i(40, 48, 52));
    for(int i = 0; i < 512; ++i) {
      auto data = octant.data(i);
      ASSERT_EQ(data.x, 5.f);
      ASSERT_EQ(data.y, 2.);
    }
  }
}

TEST(SerialiseUnitTest, SerialiseTree) {
  Octree<testT> tree;
  tree.init(1024, 10);
  const int side = se::VoxelBlock<testT>::side;
  const int max_depth = log2(tree.size());
  const int leaves_level = max_depth - log2(side);
  std::mt19937 gen(1); //Standard mersenne_twister_engine seeded with constant
  std::uniform_int_distribution<> dis(0, 1023);
  
  int num_tested = 0;
  for(int i = 1, edge = tree.size()/2; i <= leaves_level; ++i, edge = edge/2) { 
    for(int j = 0; j < 20; ++j) {
      Eigen::Vector3i vox(dis(gen), dis(gen), dis(gen));
      tree.insert(vox(0), vox(1), vox(2), i);
    }
  }
  std::string filename = "octree-test.bin";
  tree.save(filename);

  Octree<testT> tree_copy;
  tree_copy.load(filename);

  auto& node_buffer_base = tree.getNodesBuffer();
  auto& node_buffer_copy = tree_copy.getNodesBuffer();
  ASSERT_EQ(node_buffer_base.size(), node_buffer_copy.size());
  for(int i = 0; i < node_buffer_base.size(); ++i) {
    se::Node<testT> * n  = node_buffer_base[i];
    se::Node<testT> * n1 = node_buffer_copy[i];
    ASSERT_EQ(n->code, n1->code);
    ASSERT_EQ(n->children_mask_, n1->children_mask_);
  }

  auto& block_buffer_base = tree.getBlockBuffer();
  auto& block_buffer_copy = tree_copy.getBlockBuffer();
  ASSERT_EQ(block_buffer_base.size(), block_buffer_copy.size());
}

