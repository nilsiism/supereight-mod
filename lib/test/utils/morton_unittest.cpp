#include <random>
#include "math_utils.h"
#include "utils/morton_utils.hpp"
#include "octree_defines.h"
#include "gtest/gtest.h"

TEST(MortonCoding, RandomInts) {

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dis(0, 4096);

  for(int i = 0; i < 1000; ++i) {
    const Eigen::Vector3i vox = {dis(gen), dis(gen), dis(gen)};
    const octlib::key_t code = compute_morton(vox(0), vox(1), vox(2));
    const Eigen::Vector3i decoded = unpack_morton(code);
    ASSERT_EQ(decoded(0), vox(0));
    ASSERT_EQ(decoded(1), vox(1));
    ASSERT_EQ(decoded(2), vox(2));
  }

}

TEST(MortonCoding, ExhaustiveTest) {

  for(int z = 2048; z < 4096; ++z)
    for(int y = 2048; y < 2050; ++y)
      for(int x = 0; x < 4096; ++x){
        const Eigen::Vector3i vox = {x, y, z};
        const octlib::key_t code = compute_morton(vox(0), vox(1), vox(2));
        const Eigen::Vector3i decoded = unpack_morton(code);
        ASSERT_EQ(decoded(0), vox(0));
        ASSERT_EQ(decoded(1), vox(1));
        ASSERT_EQ(decoded(2), vox(2));
  }

}
