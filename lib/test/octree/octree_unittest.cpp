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
  octlib::key_t code = compute_morton(octant.x, octant.y, octant.z) | 5;
  octlib::key_t ancestor = compute_morton(64, 64, 128) | 3;
  ASSERT_EQ(true , descendant(code, ancestor, max_depth)); 

  ancestor = compute_morton(128, 64, 64) | 3;
  ASSERT_FALSE(descendant(code, ancestor, max_depth)); 
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

TEST(Octree, FarCorner) {
  /*
   * The far corner should always be "exterior", meaning that moving one step
   * in the outward direction (i.e. away from the center) in *any* direction 
   * implies leaving the parent octant. For simplicity here the corners
   * individually, but this can be done programmatically testing this property.
   * TODO: change this test case to be more exhaustive.
   */

  const int max_depth = 5;
  const int level = 2;

  /* First child */
  const octlib::key_t cell0 = compute_morton(16, 16, 16);
  const int3 fc0 = far_corner(cell0, level, max_depth);
  ASSERT_EQ(fc0.x, 16);
  ASSERT_EQ(fc0.y, 16);
  ASSERT_EQ(fc0.z, 16);

  /* Second child */
  const octlib::key_t cell1 = compute_morton(24, 16, 16);
  const int3 fc1 = far_corner(cell1, level, max_depth);
  ASSERT_EQ(fc1.x, 32);
  ASSERT_EQ(fc1.y, 16);
  ASSERT_EQ(fc1.z, 16);

  /* Third child */
  const octlib::key_t cell2 = compute_morton(16, 24, 16);
  const int3 fc2 = far_corner(cell2, level, max_depth);
  ASSERT_EQ(fc2.x, 16);
  ASSERT_EQ(fc2.y, 32);
  ASSERT_EQ(fc2.z, 16);

  /* Fourth child */
  const octlib::key_t cell3 = compute_morton(24, 24, 16);
  const int3 fc3 = far_corner(cell3, level, max_depth);
  ASSERT_EQ(fc3.x, 32);
  ASSERT_EQ(fc3.y, 32);
  ASSERT_EQ(fc3.z, 16);

  /* Fifth child */
  const octlib::key_t cell4 = compute_morton(24, 24, 16);
  const int3 fc4 = far_corner(cell4, level, max_depth);
  ASSERT_EQ(fc4.x, 32);
  ASSERT_EQ(fc4.y, 32);
  ASSERT_EQ(fc4.z, 16);

  /* sixth child */
  const octlib::key_t cell5 = compute_morton(16, 16, 24);
  const int3 fc5 = far_corner(cell5, level, max_depth);
  ASSERT_EQ(fc5.x, 16);
  ASSERT_EQ(fc5.y, 16);
  ASSERT_EQ(fc5.z, 32);

  /* seventh child */
  const octlib::key_t cell6 = compute_morton(24, 16, 24);
  const int3 fc6 = far_corner(cell6, level, max_depth);
  ASSERT_EQ(fc6.x, 32);
  ASSERT_EQ(fc6.y, 16);
  ASSERT_EQ(fc6.z, 32);

  /* eight child */
  const octlib::key_t cell7 = compute_morton(24, 24, 24);
  const int3 fc7 = far_corner(cell7, level, max_depth);
  ASSERT_EQ(fc7.x, 32);
  ASSERT_EQ(fc7.y, 32);
  ASSERT_EQ(fc7.z, 32);
}

TEST(Octree, InnerOctantExteriorNeighbours) {
  const int max_depth = 5;
  const int level = 2;
  const int side = 1 << (max_depth - level);
  const octlib::key_t cell = (compute_morton(16, 16, 16) & ~SCALE_MASK) | level;
  octlib::key_t N[7];
  exterior_neighbours(N, cell, level, max_depth);
  const octlib::key_t p = parent(cell, max_depth);
  
  const octlib::key_t neighbours_gt[7] = 
    {(compute_morton(15, 16, 16) & ~SCALE_MASK) | level,
     (compute_morton(16, 15, 16) & ~SCALE_MASK) | level,
     (compute_morton(15, 15, 16)  & ~SCALE_MASK) | level,
     (compute_morton(16, 16, 15) & ~SCALE_MASK) | level,
     (compute_morton(15, 16, 15)  & ~SCALE_MASK) | level,
     (compute_morton(16, 15, 15)  & ~SCALE_MASK) | level,
     (compute_morton(15, 15, 15)   & ~SCALE_MASK) | level };
  for(int i = 0; i < 7; ++i) {
    // std::bitset<64> c(N[i]);
    // std::bitset<64> a(p);
    // std::cout << a << std::endl;
    // std::cout << c << std::endl << std::endl;
    ASSERT_EQ(neighbours_gt[i], N[i]);
    ASSERT_FALSE(parent(N[i], max_depth) == p);
    // std::cout << (unpack_morton(N[i] & ~SCALE_MASK)) << std::endl;
  }
}

TEST(Octree, EdgeOctantExteriorNeighbours) {
  const int max_depth = 5;
  const uint size = std::pow(2, 5);
  const int level = 2;
  const octlib::key_t cell = (compute_morton(0, 16, 16) & ~SCALE_MASK) | level;
  octlib::key_t N[7];
  exterior_neighbours(N, cell, level, max_depth);
  const octlib::key_t p = parent(cell, max_depth);
  
  for(int i = 0; i < 7; ++i) {
    const uint3 corner = unpack_morton(N[i] & ~SCALE_MASK);
    ASSERT_TRUE(in(corner.x, 0u, size - 1) && 
                in(corner.y, 0u, size - 1) &&
                in(corner.z, 0u, size - 1));
  }
}
