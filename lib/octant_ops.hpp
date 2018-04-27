#ifndef OCTANT_OPS_HPP
#define OCTANT_OPS_HPP
#include "utils/morton_utils.hpp"
#include "octree_defines.h"
#include <iostream>
#include <bitset>

/*
 * Algorithm 5 of p4est paper: https://epubs.siam.org/doi/abs/10.1137/100791634
 */
inline uint3 face_neighbour(const octlib::key_t o, 
    const unsigned int face, const unsigned int l, 
    const unsigned int max_depth) {
  uint3 coords = unpack_morton(o & ~SCALE_MASK);
  const unsigned int side = 1 << (max_depth - l); 
  coords.x = coords.x + ((face == 0) ? -side : (face == 1) ? side : 0);
  coords.y = coords.y + ((face == 2) ? -side : (face == 3) ? side : 0);
  coords.z = coords.z + ((face == 4) ? -side : (face == 5) ? side : 0);
  return {coords.x, coords.y, coords.z};
}

/*
 * \brief Return true if octant is a descendant of ancestor
 * \param octant 
 * \param ancestor 
 */
inline bool descendant(octlib::key_t octant, 
    octlib::key_t ancestor, const int max_depth) {
  const int level = (ancestor & SCALE_MASK) - 1;
  const int idx = level + MAX_BITS - 1 - max_depth;
  octant = octant & MASK[idx];
  ancestor = ancestor & MASK[idx];
  return (ancestor ^ octant) == 0;
}

/*
 * \brief Computes the parent's morton code of a given octant
 * \param octant
 * \param max_depth max depth of the tree on which the octant lives
 */
inline octlib::key_t parent(const octlib::key_t& octant, const int max_depth) {
  const int level = (octant & SCALE_MASK) - 1;
  const int idx = level + MAX_BITS - 1 - max_depth;
  return (octant & MASK[idx]) | level;
}

/*
 * \brief Computes the octants's id in its local brotherhood
 * \param octant
 * \param level of octant 
 * \param max_depth max depth of the tree on which the octant lives
 */
inline int child_id(octlib::key_t octant, const int level, 
    const int max_depth) {
  int shift = max_depth - level;
  octant = (octant & ~SCALE_MASK) >> shift*3;
  int idx = (octant & 0x01) | (octant & 0x02) | (octant & 0x04);
  return idx;
}

/*
 * \brief Computes the octants's corner which is not shared with its siblings
 * \param octant
 * \param level of octant 
 * \param max_depth max depth of the tree on which the octant lives
 */
inline int3 far_corner(const octlib::key_t octant, const int level, 
    const int max_depth) {
  const unsigned int side = 1 << (max_depth - level); 
  const int idx = child_id(octant, level, max_depth);
  const uint3 coordinates = unpack_morton(octant & ~SCALE_MASK);
  return make_int3(coordinates.x + (idx & 1) * side,
                   coordinates.y + ((idx & 2) >> 1) * side,
                   coordinates.z + ((idx & 4) >> 2) * side);
}

/*
 * \brief Computes the non-sibling neighbourhood around an octants. In the
 * special case in which the octant lies on an edge, neighbour are duplicated 
 * as movement outside the enclosing cube is forbidden.
 * \param result 7-vector containing the neighbours
 * \param octant
 * \param level of octant 
 * \param max_depth max depth of the tree on which the octant lives
 */
inline void exterior_neighbours(octlib::key_t result[7], 
    const octlib::key_t octant, const int level, const int max_depth) {

  const int idx = child_id(octant, level, max_depth);
  int3 dir = make_int3((idx & 1) ? 1 : -1,
                       (idx & 2) ? 1 : -1,
                       (idx & 4) ? 1 : -1);
  int3 base = far_corner(octant, level, max_depth);
  dir.x = in(base.x + dir.x , 0, std::pow(2, max_depth) - 1) ? dir.x : 0;
  dir.y = in(base.y + dir.y , 0, std::pow(2, max_depth) - 1) ? dir.y : 0;
  dir.z = in(base.z + dir.z , 0, std::pow(2, max_depth) - 1) ? dir.z : 0;

  result[0] = (compute_morton(base.x + dir.x, base.y +     0, base.z +      0)
               & ~SCALE_MASK) | level;
  result[1] = (compute_morton(base.x +     0, base.y + dir.y, base.z +      0)
               & ~SCALE_MASK) | level; 
  result[2] = (compute_morton(base.x + dir.x, base.y + dir.y, base.z +      0)
               & ~SCALE_MASK) | level; 
  result[3] = (compute_morton(base.x +     0, base.y +     0, base.z +  dir.z)
               & ~SCALE_MASK) | level; 
  result[4] = (compute_morton(base.x + dir.x, base.y +     0, base.z +  dir.z)
               & ~SCALE_MASK) | level; 
  result[5] = (compute_morton(base.x +     0, base.y + dir.y, base.z +  dir.z)
               & ~SCALE_MASK) | level; 
  result[6] = (compute_morton(base.x + dir.x, base.y + dir.y, base.z +  dir.z)
               & ~SCALE_MASK) | level; 
}

/*
 * \brief Computes the morton number of all siblings around an octant,
 * including itself.
 * \param result 8-vector containing the neighbours
 * \param octant
 * \param max_depth max depth of the tree on which the octant lives
 */
inline void siblings(octlib::key_t result[8], 
    const octlib::key_t octant, const int max_depth) {
  const int level = (octant & SCALE_MASK);
  const int shift = 3*(max_depth - level);
  const octlib::key_t p = parent(octant, max_depth) + 1; // set-up next level
  for(int i = 0; i < 8; ++i) {
    result[i] = p | (i << shift);
  }
}
#endif

