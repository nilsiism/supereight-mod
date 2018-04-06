#include <random>
#include "math_utils.h"
#include "utils/morton_utils.hpp"
#include "octree_defines.h"
#include "gtest/gtest.h"

TEST(MortonCoding, RandomInts) {

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<unsigned int> dis(0, 4096);

  for(int i = 0; i < 1000; ++i) {
    const uint3 vox = {dis(gen), dis(gen), dis(gen)};
    const octlib::key_t code = compute_morton(vox.x, vox.y, vox.z);
    const uint3 decoded = unpack_morton(code);
    ASSERT_EQ(decoded.x, vox.x);
    ASSERT_EQ(decoded.y, vox.y);
    ASSERT_EQ(decoded.z, vox.z);
  }

}

TEST(MortonCoding, ExhaustiveTest) {

  for(unsigned int z = 2048; z < 4096; ++z)
    for(unsigned int y = 2048; y < 2050; ++y)
      for(unsigned int x = 0; x < 4096; ++x){
        const uint3 vox = {x, y, z};
        const octlib::key_t code = compute_morton(vox.x, vox.y, vox.z);
        const uint3 decoded = unpack_morton(code);
        ASSERT_EQ(decoded.x, vox.x);
        ASSERT_EQ(decoded.y, vox.y);
        ASSERT_EQ(decoded.z, vox.z);
  }
}

TEST(MortonCoding, FaceNeighbours) {
  const uint3 octant = {112, 80, 160};
  const octlib::key_t code = compute_morton(octant.x, octant.y, octant.z);
  const unsigned int max_depth = 8;
  const unsigned int leaves_depth = 5;
  const unsigned int side = 8;
  const int3 faces[6] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, 
    {0, 0, -1}, {0, 0, 1}};
  for(int i = 0; i < 6; ++i) {
    const int3 neighbour = make_int3(octant) + side * faces[i];
    const uint3 computed = face_neighbour(code, i, leaves_depth, max_depth); 
    ASSERT_EQ(neighbour.x, computed.x); 
    ASSERT_EQ(neighbour.y, computed.y); 
    ASSERT_EQ(neighbour.z, computed.z); 
  }
}
