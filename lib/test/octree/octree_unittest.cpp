#include "utils/eigen_helper.h"
#include "gtest/gtest.h"
#include "octant_ops.hpp"
#include <bitset>

TEST(Octree, OctantFaceNeighbours) {
  const Eigen::Vector3i octant = {112, 80, 160};
  const unsigned int max_depth = 8;
  const unsigned int leaves_depth = 5;
  const octlib::key_t code = 
    octlib::keyops::encode(octant(0), octant(1), octant(2), leaves_depth, max_depth);
  const unsigned int side = 8;
  const Eigen::Vector3i faces[6] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, 
    {0, 0, -1}, {0, 0, 1}};
  for(int i = 0; i < 6; ++i) {
    const Eigen::Vector3i neighbour = octant + side * faces[i];
    const Eigen::Vector3i computed = face_neighbour(code, i, leaves_depth, max_depth); 
    ASSERT_EQ(neighbour(0), computed(0)); 
    ASSERT_EQ(neighbour(1), computed(1)); 
    ASSERT_EQ(neighbour(2), computed(2)); 
  }
}

TEST(Octree, OctantDescendant) {
  const unsigned max_depth = 8;
  Eigen::Vector3i octant = {110, 80, 159};
  octlib::key_t code = 
    octlib::keyops::encode(octant(0), octant(1), octant(2), 5, max_depth);
  octlib::key_t ancestor = 
    octlib::keyops::encode(96, 64, 128, 3, max_depth);
  ASSERT_EQ(true , descendant(code, ancestor, max_depth)); 

  ancestor = octlib::keyops::encode(128, 64, 64, 3, max_depth);
  ASSERT_FALSE(descendant(code, ancestor, max_depth)); 
}

TEST(Octree, OctantParent) {
  const int max_depth = 8;
  Eigen::Vector3i octant = {112, 80, 160};
  octlib::key_t code = 
    octlib::keyops::encode(octant(0), octant(1), octant(2), 5, max_depth);
  octlib::key_t p = parent(code, max_depth);
  ASSERT_EQ(octlib::keyops::code(code), octlib::keyops::code(p));
  ASSERT_EQ(4, p & SCALE_MASK);

  code = p;
  p = parent(code, max_depth); 
  ASSERT_EQ(3, octlib::keyops::level(p));
  ASSERT_EQ(p, octlib::keyops::encode(96, 64, 160, 3, max_depth));

  code = p;
  p = parent(code, max_depth); 
  ASSERT_EQ(2, octlib::keyops::level(p));
  ASSERT_EQ(p, octlib::keyops::encode(64, 64, 128, 2, max_depth));
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
  const octlib::key_t cell0 = 
    octlib::keyops::encode(16, 16, 16, level, max_depth);
  const Eigen::Vector3i fc0 = far_corner(cell0, level, max_depth);
  ASSERT_EQ(fc0(0), 16);
  ASSERT_EQ(fc0(1), 16);
  ASSERT_EQ(fc0(2), 16);

  /* Second child */
  const octlib::key_t cell1 = octlib::keyops::encode(24, 16, 16, level, max_depth);
  const Eigen::Vector3i fc1 = far_corner(cell1, level, max_depth);
  ASSERT_EQ(fc1(0), 32);
  ASSERT_EQ(fc1(1), 16);
  ASSERT_EQ(fc1(2), 16);

  /* Third child */
  const octlib::key_t cell2 = octlib::keyops::encode(16, 24, 16, level, max_depth);
  const Eigen::Vector3i fc2 = far_corner(cell2, level, max_depth);
  ASSERT_EQ(fc2(0), 16);
  ASSERT_EQ(fc2(1), 32);
  ASSERT_EQ(fc2(2), 16);

  /* Fourth child */
  const octlib::key_t cell3 = octlib::keyops::encode(24, 24, 16, level, max_depth);
  const Eigen::Vector3i fc3 = far_corner(cell3, level, max_depth);
  ASSERT_EQ(fc3(0), 32);
  ASSERT_EQ(fc3(1), 32);
  ASSERT_EQ(fc3(2), 16);

  /* Fifth child */
  const octlib::key_t cell4 = octlib::keyops::encode(24, 24, 16, level, max_depth);
  const Eigen::Vector3i fc4 = far_corner(cell4, level, max_depth);
  ASSERT_EQ(fc4(0), 32);
  ASSERT_EQ(fc4(1), 32);
  ASSERT_EQ(fc4(2), 16);

  /* sixth child */
  const octlib::key_t cell5 = octlib::keyops::encode(16, 16, 24, level, max_depth);
  const Eigen::Vector3i fc5 = far_corner(cell5, level, max_depth);
  ASSERT_EQ(fc5(0), 16);
  ASSERT_EQ(fc5(1), 16);
  ASSERT_EQ(fc5(2), 32);

  /* seventh child */
  const octlib::key_t cell6 = octlib::keyops::encode(24, 16, 24, level, max_depth);
  const Eigen::Vector3i fc6 = far_corner(cell6, level, max_depth);
  ASSERT_EQ(fc6(0), 32);
  ASSERT_EQ(fc6(1), 16);
  ASSERT_EQ(fc6(2), 32);

  /* eight child */
  const octlib::key_t cell7 = octlib::keyops::encode(24, 24, 24, level, max_depth);
  const Eigen::Vector3i fc7 = far_corner(cell7, level, max_depth);
  ASSERT_EQ(fc7(0), 32);
  ASSERT_EQ(fc7(1), 32);
  ASSERT_EQ(fc7(2), 32);
}

TEST(Octree, InnerOctantExteriorNeighbours) {
  const int max_depth = 5;
  const int level = 2;
  const int side = 1 << (max_depth - level);
  const octlib::key_t cell = octlib::keyops::encode(16, 16, 16, level, max_depth);
  octlib::key_t N[7];
  exterior_neighbours(N, cell, level, max_depth);
  const octlib::key_t p = parent(cell, max_depth);
  
  const octlib::key_t neighbours_gt[7] = 
    {octlib::keyops::encode(15, 16, 16, level, max_depth),
     octlib::keyops::encode(16, 15, 16, level, max_depth),
     octlib::keyops::encode(15, 15, 16, level, max_depth),
     octlib::keyops::encode(16, 16, 15, level, max_depth),
     octlib::keyops::encode(15, 16, 15, level, max_depth),
     octlib::keyops::encode(16, 15, 15, level, max_depth),
     octlib::keyops::encode(15, 15, 15, level, max_depth)};
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
  const int size = 1 << max_depth;
  const int level = 2;
  const octlib::key_t cell = octlib::keyops::encode(0, 16, 16, level, max_depth);
  octlib::key_t N[7];
  exterior_neighbours(N, cell, level, max_depth);
  const octlib::key_t p = parent(cell, max_depth);
  
  for(int i = 0; i < 7; ++i) {
    const Eigen::Vector3i corner = unpack_morton(N[i] & ~SCALE_MASK);
    const int res = ((corner.array() >= Eigen::Vector3i::Constant(0).array()) 
     && (corner.array() <= Eigen::Vector3i::Constant(size - 1).array())).all();
    ASSERT_TRUE(res);
  }
}

TEST(Octree, OctantSiblings) {
  const int max_depth = 5;
  const uint size = std::pow(2, 5);
  const int level = 2;
  const octlib::key_t cell = octlib::keyops::encode(16, 16, 16, level, max_depth);
  octlib::key_t s[8];
  siblings(s, cell, max_depth);

  const int childidx = child_id(cell, level, max_depth);
  ASSERT_EQ(s[childidx], cell);
  
  for(int i = 0; i < 8; ++i) {
    // std::cout << (unpack_morton(s[i] & ~SCALE_MASK)) << std::endl;
    ASSERT_TRUE(parent(s[i], max_depth) == parent(cell, max_depth));
  }
}
