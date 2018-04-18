#include "math_utils.h"
#include "gtest/gtest.h"
#include "octant_ops.hpp"

TEST(Octree, OctantFaceNeighbours) {
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

TEST(Octree, OctantDescendant) {
  const unsigned max_depth = 8;
  uint3 octant = {112, 80, 160};
  octlib::key_t code = compute_morton(octant.x, octant.y, octant.z);
  octlib::key_t root = 0;
  ASSERT_EQ(root, descendant(code, root)); 
}
