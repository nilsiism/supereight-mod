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
  octlib::key_t ancestor = compute_morton(64, 64, 128);
  ASSERT_EQ(true , descendant(code, ancestor)); 
}

TEST(Octree, OctantParent) {
  const int max_depth = 8;
  uint3 octant = {112, 80, 160};
  octlib::key_t code = compute_morton(octant.x, octant.y, octant.z) | 5;
  octlib::key_t p = parent(code, max_depth);
  ASSERT_EQ(code & ~SCALE_MASK, p & ~SCALE_MASK);
  ASSERT_EQ(4, p & SCALE_MASK);

  code = p;
  p = parent(code, max_depth); 
  ASSERT_EQ(3, p & SCALE_MASK);
  ASSERT_EQ(p & ~SCALE_MASK, compute_morton(96, 64, 160));

  code = p;
  p = parent(code, max_depth); 
  ASSERT_EQ(2, p & SCALE_MASK);
  ASSERT_EQ(p & ~SCALE_MASK, compute_morton(64, 64, 128));
}
