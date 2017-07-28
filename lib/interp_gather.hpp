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

template <typename FieldType, template<typename FieldType> class MapIndex,
         class FieldSelector>
inline void gather_points(const MapIndex<FieldType>& fetcher, const uint3 base, 
    FieldSelector select, float points[8]) {
 
  unsigned int blockSize =  VoxelBlock<FieldType>::side;
  unsigned int crossmask = ((base.x % blockSize) == blockSize - 1 << 2) | 
                           ((base.y % blockSize) == blockSize - 1 << 1) |
                            (base.z % blockSize) == blockSize - 1;

  switch(crossmask) {
    case 0: /* all local */
      {
        VoxelBlock<FieldType> * block = fetcher.fetch(base.x, base.y, base.z);
        gather_local(block, base, select, points);
      }
      break;
    case 1: /* z crosses */
      {
        const unsigned int offs1[4] = {0, 1, 2, 3};
        const unsigned int offs2[4] = {4, 5, 6, 7};
        VoxelBlock<FieldType> * block = fetcher.fetch(base.x, base.y, base.z);
        if(!block) {
          std::cerr << "Fetched wrong block" << std::endl;
          points[offs1[0]] = select(VoxelBlock<FieldType>::empty());
          points[offs1[1]] = select(VoxelBlock<FieldType>::empty());
          points[offs1[2]] = select(VoxelBlock<FieldType>::empty());
          points[offs1[3]] = select(VoxelBlock<FieldType>::empty());
        }
        gather_4(block, base, [](const auto& val){ return val; }, offs1, points);
        const uint3 base1 = base + offs2[0];
        block = fetcher.fetch(base1.x, base1.y, base1.z);
        if(!block) {
          std::cerr << "Fetched wrong block" << std::endl;
          points[offs2[0]] = select(VoxelBlock<FieldType>::empty());
          points[offs2[1]] = select(VoxelBlock<FieldType>::empty());
          points[offs2[2]] = select(VoxelBlock<FieldType>::empty());
          points[offs2[3]] = select(VoxelBlock<FieldType>::empty());
        }
        gather_4(block, base, [](const auto& val){ return val; }, offs2, points);
      }
      break;
    case 2: /* y crosses */ 
      {
        const unsigned int offs1[4] = {0, 1, 4, 5};
        const unsigned int offs2[4] = {2, 3, 6, 7};
      }
      break;
    case 3: /* y, z cross */ 
      {
        const unsigned int offs1[2] = {0, 1};
        const unsigned int offs2[2] = {2, 3};
        const unsigned int offs3[2] = {4, 5};
        const unsigned int offs4[2] = {6, 7};
      }
      break;
    case 4: /* x crosses */ 
      {
        const unsigned int offs1[4] = {0, 2, 4, 6};
        const unsigned int offs2[4] = {1, 3, 5, 7};
      }
    case 5: /* x,z cross */ 
      {
        const unsigned int offs1[2] = {0, 2};
        const unsigned int offs2[2] = {1, 3};
        const unsigned int offs3[2] = {4, 6};
        const unsigned int offs4[2] = {5, 7};
      }
// static constexpr const uint3 interp_offsets[8] = 
//   {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, 
//    {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
    case 6: /* x,y cross */ 
      {
        const unsigned int offs1[2] = {0, 4};
        const unsigned int offs2[2] = {1, 5};
        const unsigned int offs3[2] = {2, 6};
        const unsigned int offs4[2] = {3, 7};
      }
  }
}
#endif
