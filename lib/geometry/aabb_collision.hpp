#ifndef AABB_COLLISION_HPP
#define AABB_COLLISION_HPP
#include <cmath>
#include <math_utils.h> 

enum class collision_status {
  occupied,
  unseen,
  empty
};

namespace geometry {
  inline int axis_collision(int a, const int a_edge, 
       int b, const int b_edge) {
    /* Half plane intersection test */
    a = a + (a_edge/2);
    b = b + (b_edge/2);
    return (std::abs(b - a) > (a_edge + b_edge)/2) ? 0 : 1;
  }

  inline int axis_collision(float a, const float a_edge, 
       float b, const float b_edge) {
    /* Half plane intersection test */
    a = a + (a_edge/2);
    b = b + (b_edge/2);
    return (std::fabs(b - a) > (a_edge + b_edge)/2) ? 0 : 1;
  }

  inline int aabb_aabb_collision(const int3 a, const int3 a_edge, 
      const int3 b, const int3 b_edge){

    return axis_collision(a.x, a_edge.x, b.x, b_edge.x) && 
           axis_collision(a.y, a_edge.y, b.y, b_edge.y) && 
           axis_collision(a.z, a_edge.z, b.z, b_edge.z);
  }
}
#endif
