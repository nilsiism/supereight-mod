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
    const morton_type code = compute_morton(vox.x, vox.y, vox.z);
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
        const morton_type code = compute_morton(vox.x, vox.y, vox.z);
        const uint3 decoded = unpack_morton(code);
        ASSERT_EQ(decoded.x, vox.x);
        ASSERT_EQ(decoded.y, vox.y);
        ASSERT_EQ(decoded.z, vox.z);
  }

}
