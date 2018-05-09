#include "math_utils.h"
#include "geometry/aabb_collision.hpp"
#include "gtest/gtest.h"

TEST(AABBAABBTest, SquareOverlap) {

  const Eigen::Vector3i a = Eigen::Vector3i(5, 6, 9);
  const Eigen::Vector3i a_edge = Eigen::Vector3i::Constant(3);
  const Eigen::Vector3i b = Eigen::Vector3i(4, 4, 8);
  const Eigen::Vector3i b_edge = Eigen::Vector3i::Constant(3);
  int overlaps = geometry::aabb_aabb_collision(a, a_edge, b, b_edge);
  ASSERT_EQ(overlaps, 1);
}

TEST(AABBAABBTest, SquareDisjoint1Axis) {

  const Eigen::Vector3i a = Eigen::Vector3i(5, 6, 9);
  const Eigen::Vector3i a_edge = Eigen::Vector3i::Constant(3);
  const Eigen::Vector3i b = Eigen::Vector3i(4, 4, 1);
  const Eigen::Vector3i b_edge = Eigen::Vector3i::Constant(3);
  int overlaps = geometry::aabb_aabb_collision(a, a_edge, b, b_edge);
  ASSERT_EQ(overlaps, 0);
}

TEST(AABBAABBTest, SquareDisjoint2Axis) {
  /* Disjoint on y and z */
  const Eigen::Vector3i a = Eigen::Vector3i(5, 6, 9);
  const Eigen::Vector3i a_edge = Eigen::Vector3i::Constant(3);
  const Eigen::Vector3i b = Eigen::Vector3i(6, 22, 13);
  const Eigen::Vector3i b_edge = Eigen::Vector3i::Constant(10);

  int overlapx = geometry::axis_overlap(a(0), a_edge(0), b(0), b_edge(0));
  int overlapy = geometry::axis_overlap(a(1), a_edge(1), b(1), b_edge(1));
  int overlapz = geometry::axis_overlap(a(2), a_edge(2), b(2), b_edge(2));

  ASSERT_EQ(overlapx, 1);
  ASSERT_EQ(overlapy, 0);
  ASSERT_EQ(overlapz, 0);

  int overlaps = geometry::aabb_aabb_collision(a, a_edge, b, b_edge);
  ASSERT_EQ(overlaps, 0);
}

TEST(AABBAABBTest, SquareDisjoint) {
  /* Disjoint on x, y and z */
  const Eigen::Vector3i a = Eigen::Vector3i(5, 6, 9);
  const Eigen::Vector3i a_edge = Eigen::Vector3i::Constant(4);
  const Eigen::Vector3i b = Eigen::Vector3i(12, 22, 43);
  const Eigen::Vector3i b_edge = Eigen::Vector3i::Constant(10);

  int overlapx = geometry::axis_overlap(a(0), a_edge(0), b(0), b_edge(0));
  int overlapy = geometry::axis_overlap(a(1), a_edge(1), b(1), b_edge(1));
  int overlapz = geometry::axis_overlap(a(2), a_edge(2), b(2), b_edge(2));

  ASSERT_EQ(overlapx, 0);
  ASSERT_EQ(overlapy, 0);
  ASSERT_EQ(overlapz, 0);

  int overlaps = geometry::aabb_aabb_collision(a, a_edge, b, b_edge);
  ASSERT_EQ(overlaps, 0);
}

TEST(AABBAABBTest, SquareEnclosed) {
  /* Disjoint on x, y and z */
  const Eigen::Vector3i a = Eigen::Vector3i(5, 6, 9);
  const Eigen::Vector3i a_edge = Eigen::Vector3i::Constant(10);
  const Eigen::Vector3i b = Eigen::Vector3i(6, 7, 8);
  const Eigen::Vector3i b_edge = Eigen::Vector3i::Constant(2);

  int overlapx = geometry::axis_overlap(a(0), a_edge(0), b(0), b_edge(0));
  int overlapy = geometry::axis_overlap(a(1), a_edge(1), b(1), b_edge(1));
  int overlapz = geometry::axis_overlap(a(2), a_edge(2), b(2), b_edge(2));

  ASSERT_EQ(overlapx, 1);
  ASSERT_EQ(overlapy, 1);
  ASSERT_EQ(overlapz, 1);

  int overlaps = geometry::aabb_aabb_collision(a, a_edge, b, b_edge);
  ASSERT_EQ(overlaps, 1);
}

TEST(AABBAABBTest, Inclusion) {
  const Eigen::Vector3i a = Eigen::Vector3i(2, 1, 3); 
  const Eigen::Vector3i a_edge = Eigen::Vector3i::Constant(10);
  const Eigen::Vector3i b = Eigen::Vector3i(3, 4, 5);
  const Eigen::Vector3i b_edge = Eigen::Vector3i::Constant(2);
  int included = geometry::aabb_aabb_inclusion(a, a_edge, b, b_edge);
  ASSERT_EQ(included, 1);
}
