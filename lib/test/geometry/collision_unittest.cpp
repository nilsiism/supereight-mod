#include "math_utils.h"
#include "geometry/aabb_collision.hpp"
#include "gtest/gtest.h"

TEST(AABBAABBTest, SquareOverlap) {

  const int3 a = make_int3(5, 6, 9);
  const int3 a_edge = make_int3(3);
  const int3 b = make_int3(4, 4, 8);
  const int3 b_edge = make_int3(3);
  int overlaps = geometry::aabb_aabb_collision(a, a_edge, b, b_edge);
  ASSERT_EQ(overlaps, 1);
}

TEST(AABBAABBTest, SquareDisjoint1Axis) {

  const int3 a = make_int3(5, 6, 9);
  const int3 a_edge = make_int3(3);
  const int3 b = make_int3(4, 4, 1);
  const int3 b_edge = make_int3(3);
  int overlaps = geometry::aabb_aabb_collision(a, a_edge, b, b_edge);
  ASSERT_EQ(overlaps, 0);
}

TEST(AABBAABBTest, SquareDisjoint2Axis) {
  /* Disjoint on y and z */
  const int3 a = make_int3(5, 6, 9);
  const int3 a_edge = make_int3(3);
  const int3 b = make_int3(6, 22, 13);
  const int3 b_edge = make_int3(10);

  int overlapx = geometry::axis_collision(a.x, a_edge.x, b.x, b_edge.x);
  int overlapy = geometry::axis_collision(a.y, a_edge.y, b.y, b_edge.y);
  int overlapz = geometry::axis_collision(a.z, a_edge.z, b.z, b_edge.z);

  ASSERT_EQ(overlapx, 1);
  ASSERT_EQ(overlapy, 0);
  ASSERT_EQ(overlapz, 0);

  int overlaps = geometry::aabb_aabb_collision(a, a_edge, b, b_edge);
  ASSERT_EQ(overlaps, 0);
}

TEST(AABBAABBTest, SquareDisjoint) {
  /* Disjoint on x, y and z */
  const int3 a = make_int3(5, 6, 9);
  const int3 a_edge = make_int3(4);
  const int3 b = make_int3(12, 22, 43);
  const int3 b_edge = make_int3(10);

  int overlapx = geometry::axis_collision(a.x, a_edge.x, b.x, b_edge.x);
  int overlapy = geometry::axis_collision(a.y, a_edge.y, b.y, b_edge.y);
  int overlapz = geometry::axis_collision(a.z, a_edge.z, b.z, b_edge.z);

  ASSERT_EQ(overlapx, 0);
  ASSERT_EQ(overlapy, 0);
  ASSERT_EQ(overlapz, 0);

  int overlaps = geometry::aabb_aabb_collision(a, a_edge, b, b_edge);
  ASSERT_EQ(overlaps, 0);
}
