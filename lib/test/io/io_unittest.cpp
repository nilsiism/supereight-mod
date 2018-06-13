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


TEST(SerialiseUnitTest, WriteReadNode) {
  std::string filename = "test.bin";
  {
    std::ofstream os (filename, std::ios::binary); 
    se::Node<testT> octant;
    octant.code = 24;
    octant.side = 256;
    se::internal::serialise(os, octant);
  }

  {
    std::ifstream is(filename, std::ios::binary);
    se::Node<testT> octant;
    se::internal::deserialise(octant, is);
    std::cout << "Read: " << octant.code << ", " << octant.side << " ";
    for(int i = 0; i < 8; ++i)
      std::cout << octant.value_[i] << " ";
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
    std::cout << "Read: " << octant.code << ", " << octant.side << " "
      << ", coords:\n " << octant.coordinates();
    for(int i = 0; i < 512; ++i)
      ASSERT_EQ(octant.data(i), 5.f);
  }
}

