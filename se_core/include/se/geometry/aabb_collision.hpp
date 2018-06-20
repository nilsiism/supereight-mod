#ifndef AABB_COLLISION_HPP
#define AABB_COLLISION_HPP
#include <cmath>
#include "../utils/se_common.h"

namespace geometry {
  inline int axis_overlap(int a, const int a_edge, 
       int b, const int b_edge) {
    /* Half plane intersection test */
    a = a + (a_edge/2);
    b = b + (b_edge/2);
    return (std::abs(b - a) > (a_edge + b_edge)/2) ? 0 : 1;
  }

  inline int axis_overlap(float a, const float a_edge, 
       float b, const float b_edge) {
    /* Half plane intersection test */
    a = a + (a_edge/2);
    b = b + (b_edge/2);
    return (std::fabs(b - a) > (a_edge + b_edge)/2) ? 0 : 1;
  }

  inline int axis_contained(float a, const float a_edge, 
       float b, const float b_edge) {
    /* Segment a includes segment b */
    return (a < b) && ((a + a_edge) > (b + b_edge)); 
  }


  inline int aabb_aabb_collision(const Eigen::Vector3i a, const Eigen::Vector3i a_edge, 
      const Eigen::Vector3i b, const Eigen::Vector3i b_edge){

    return axis_overlap(a(0), a_edge(0), b(0), b_edge(0)) && 
           axis_overlap(a(1), a_edge(1), b(1), b_edge(1)) && 
           axis_overlap(a(2), a_edge(2), b(2), b_edge(2));
  }

  inline int aabb_aabb_inclusion(const Eigen::Vector3i a, const Eigen::Vector3i a_edge, 
      const Eigen::Vector3i b, const Eigen::Vector3i b_edge){
    /* Box a contains box b */
    return axis_contained(a(0), a_edge(0), b(0), b_edge(0)) && 
           axis_contained(a(1), a_edge(1), b(1), b_edge(1)) && 
           axis_contained(a(2), a_edge(2), b(2), b_edge(2));
  }
}
#endif
