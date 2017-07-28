#ifndef INTERP_GATHER_H
#define INTERP_GATHER_H
#include "node.hpp"

/*
 * Interpolation's point gather offsets
 */

static constexpr const uint3 interp_offsets[8] = 
  {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, 
   {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};

template <typename FieldType, typename FieldSelector>
inline void gather_local(const VoxelBlock<FieldType>* block, const uint3 base, 
    FieldSelector select, float points[8]) {

  points[0] = select(block->data(base + interp_offsets[0]));
  points[1] = select(block->data(base + interp_offsets[1]));
  points[2] = select(block->data(base + interp_offsets[2]));
  points[3] = select(block->data(base + interp_offsets[3]));
  points[4] = select(block->data(base + interp_offsets[4]));
  points[5] = select(block->data(base + interp_offsets[5]));
  points[6] = select(block->data(base + interp_offsets[6]));
  points[7] = select(block->data(base + interp_offsets[7]));

}

template <typename FieldType, typename FieldSelector>
inline void gather_4(const VoxelBlock<FieldType>* block, const uint3 base, 
    FieldSelector select, const unsigned int offsets[4], float points[8]) {

  points[offsets[0]] = select(block->data(base + interp_offsets[offsets[0]]));
  points[offsets[1]] = select(block->data(base + interp_offsets[offsets[1]]));
  points[offsets[2]] = select(block->data(base + interp_offsets[offsets[2]]));
  points[offsets[3]] = select(block->data(base + interp_offsets[offsets[3]]));
}

#endif
