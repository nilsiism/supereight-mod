#ifndef OCTANT_OPS_HPP
#define OCTANT_OPS_HPP
#include "utils/morton_utils.hpp"
#include "octree_defines.h"
/*
 * Algorithm 5 of p4est paper: https://epubs.siam.org/doi/abs/10.1137/100791634
 */
inline uint3 face_neighbour(const octlib::key_t o, 
    const unsigned int face, const unsigned int l, 
    const unsigned int max_depth) {
  uint3 coords = unpack_morton(o);
  const unsigned int side = 1 << (max_depth - l); 
  coords.x = coords.x + ((face == 0) ? -side : (face == 1) ? side : 0);
  coords.y = coords.y + ((face == 2) ? -side : (face == 3) ? side : 0);
  coords.z = coords.z + ((face == 4) ? -side : (face == 5) ? side : 0);
  return {coords.x, coords.y, coords.z};
}

#endif

