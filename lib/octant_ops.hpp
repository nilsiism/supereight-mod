#ifndef OCTANT_OPS_HPP
#define OCTANT_OPS_HPP
#include "utils/morton_utils.hpp"
#include "octree_defines.h"
#include <bitset>
#include <iostream>

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
inline bool descendant(const octlib::key_t& octant,
    const octlib::key_t& ancestor) {
  return (ancestor & octant) == ancestor;
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
  octant = octant >> shift*3;
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
#endif

