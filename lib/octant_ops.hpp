#ifndef OCTANT_OPS_HPP
#define OCTANT_OPS_HPP
#include "utils/morton_utils.hpp"
#include "octree_defines.h"
#include <iostream>
#include <bitset>
#include <Eigen/Dense>

namespace octlib {
  namespace keyops {

    octlib::key_t code(const octlib::key_t key) {
      return key & ~SCALE_MASK;
    }

    int level(const octlib::key_t key) {
      return key & SCALE_MASK;
    }

    octlib::key_t encode(const int x, const int y, const int z, 
        const int level, const int max_depth) {
      const int offset = MAX_BITS - max_depth + level - 1;
      return (compute_morton(x, y, z) & MASK[offset]) | level;
    }

    Eigen::Vector3i decode(const octlib::key_t key) {
      return unpack_morton(key & ~SCALE_MASK);
    }
  }
}

/*
 * Algorithm 5 of p4est paper: https://epubs.siam.org/doi/abs/10.1137/100791634
 */
inline Eigen::Vector3i face_neighbour(const octlib::key_t o, 
    const unsigned int face, const unsigned int l, 
    const unsigned int max_depth) {
  Eigen::Vector3i coords = octlib::keyops::decode(o);
  const unsigned int side = 1 << (max_depth - l); 
  coords(0) = coords(0) + ((face == 0) ? -side : (face == 1) ? side : 0);
  coords(1) = coords(1) + ((face == 2) ? -side : (face == 3) ? side : 0);
  coords(2) = coords(2) + ((face == 4) ? -side : (face == 5) ? side : 0);
  return {coords(0), coords(1), coords(2)};
}

/*
 * \brief Return true if octant is a descendant of ancestor
 * \param octant 
 * \param ancestor 
 */
inline bool descendant(octlib::key_t octant, 
    octlib::key_t ancestor, const int max_depth) {
  const int level = octlib::keyops::level(ancestor);
  const int idx = MAX_BITS - max_depth - 1 + level;
  octant = octant & MASK[idx];
  ancestor = octlib::keyops::code(ancestor);
  return (ancestor ^ octant) == 0;
}

/*
 * \brief Computes the parent's morton code of a given octant
 * \param octant
 * \param max_depth max depth of the tree on which the octant lives
 */
inline octlib::key_t parent(const octlib::key_t& octant, const int max_depth) {
  const int level = octlib::keyops::level(octant) - 1;
  const int idx = MAX_BITS - max_depth + level - 1;
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
  octant = octlib::keyops::code(octant) >> shift*3;
  int idx = (octant & 0x01) | (octant & 0x02) | (octant & 0x04);
  return idx;
}

/*
 * \brief Computes the octants's corner which is not shared with its siblings
 * \param octant
 * \param level of octant 
 * \param max_depth max depth of the tree on which the octant lives
 */
inline Eigen::Vector3i far_corner(const octlib::key_t octant, const int level, 
    const int max_depth) {
  const unsigned int side = 1 << (max_depth - level); 
  const int idx = child_id(octant, level, max_depth);
  const Eigen::Vector3i coordinates = octlib::keyops::decode(octant);
  return Eigen::Vector3i(coordinates(0) + (idx & 1) * side,
                   coordinates(1) + ((idx & 2) >> 1) * side,
                   coordinates(2) + ((idx & 4) >> 2) * side);
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
  Eigen::Vector3i dir = Eigen::Vector3i((idx & 1) ? 1 : -1,
                       (idx & 2) ? 1 : -1,
                       (idx & 4) ? 1 : -1);
  Eigen::Vector3i base = far_corner(octant, level, max_depth);
  dir(0) = in(base(0) + dir(0) , 0, std::pow(2, max_depth) - 1) ? dir(0) : 0;
  dir(1) = in(base(1) + dir(1) , 0, std::pow(2, max_depth) - 1) ? dir(1) : 0;
  dir(2) = in(base(2) + dir(2) , 0, std::pow(2, max_depth) - 1) ? dir(2) : 0;

 result[0] = octlib::keyops::encode(base(0) + dir(0), base(1) + 0, base(2) + 0, 
     level, max_depth);
 result[1] = octlib::keyops::encode(base(0) + 0, base(1) + dir(1), base(2) + 0, 
     level, max_depth); 
 result[2] = octlib::keyops::encode(base(0) + dir(0), base(1) + dir(1), base(2) + 0, 
     level, max_depth); 
 result[3] = octlib::keyops::encode(base(0) + 0, base(1) + 0, base(2) + dir(2), 
     level, max_depth); 
 result[4] = octlib::keyops::encode(base(0) + dir(0), base(1) + 0, base(2) + dir(2), 
     level, max_depth); 
 result[5] = octlib::keyops::encode(base(0) + 0, base(1) + dir(1), base(2) + dir(2), 
     level, max_depth); 
 result[6] = octlib::keyops::encode(base(0) + dir(0), base(1) + dir(1), 
     base(2) + dir(2), level, max_depth); 
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

