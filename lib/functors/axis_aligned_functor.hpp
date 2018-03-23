#ifndef AA_FUNCTOR_HPP
#define AA_FUNCTOR_HPP
#include <functional>
#include <vector>

#include "math_utils.h"
#include "algorithms/filter.hpp"
#include "node.hpp"
#include "functors/data_handler.hpp"

namespace iterators {
  namespace functor {

  template <typename FieldType, template <typename FieldT> class MapT, 
            typename UpdateF>

    class axis_aligned {
      public:
      axis_aligned(MapT<FieldType>& map, UpdateF f) : _map(map), _function(f),
      _min(make_int3(0)), _max(make_int3(map.size())){ }

      void update_block(VoxelBlock<FieldType> * block) {
        const int3 blockCoord = block->coordinates();
        unsigned int y, z, x; 
        unsigned int blockSide = VoxelBlock<FieldType>::side;
        unsigned int xlast = blockCoord.x + blockSide;
        unsigned int ylast = blockCoord.y + blockSide;
        unsigned int zlast = blockCoord.z + blockSide;

        for(z = blockCoord.z; z < zlast; ++z) {
          for (y = blockCoord.y; y < ylast; ++y) {
            for (x = blockCoord.x; x < xlast; ++x) {
              int3 vox = make_int3(x, y, z);
              VoxelBlockHandler<FieldType> handler = {block, vox};
              _function(handler, vox);
            }
          }
        }
      }

      void update_node(Node<FieldType> * node) { 
        const int3 voxel = make_int3(unpack_morton(node->code));
        for(int i = 0; i < 8; ++i) {
          const int3 dir =  make_int3((i & 1) > 0, (i & 2) > 0, (i & 4) > 0);
          NodeHandler<FieldType> handler = {node, i};
          _function(handler, voxel + (dir * (node->side/2)));
        }
      }

      void apply() {

        auto& block_list = _map.getBlockBuffer();
        size_t list_size = block_list.size();
#pragma omp parallel for
        for(unsigned int i = 0; i < list_size; ++i){
          update_block(block_list[i]);
        }

        auto& nodes_list = _map.getNodesBuffer();
        list_size = nodes_list.size();
#pragma omp parallel for
        for(unsigned int i = 0; i < list_size; ++i){
          update_node(nodes_list[i]);
        }
      }

    private:
      MapT<FieldType>& _map; 
      UpdateF _function; 
      int3 _min;
      int3 _max;
    };
  }
}
#endif
