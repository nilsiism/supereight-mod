#include <fstream>
#include <string>
#include "se_serialise.hpp"
#include "node.hpp"
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

