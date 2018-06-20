#ifndef AA_FUNCTOR_HPP
#define AA_FUNCTOR_HPP
#include <functional>
#include <vector>

#include "utils/se_common.h"
#include "algorithms/filter.hpp"
#include "node.hpp"
#include "functors/data_handler.hpp"
#include "geometry/aabb_collision.hpp"

namespace se {
  namespace functor {

    template <typename FieldType, template <typename FieldT> class MapT, 
              typename UpdateF>

      class axis_aligned {
        public:
        axis_aligned(MapT<FieldType>& map, UpdateF f) : _map(map), _function(f),
        _min(Eigen::Vector3i::Constant(0)), 
        _max(Eigen::Vector3i::Constant(map.size())){ }

        axis_aligned(MapT<FieldType>& map, UpdateF f, const Eigen::Vector3i min,
            const Eigen::Vector3i max) : _map(map), _function(f),
        _min(min), _max(max){ }

        void update_block(se::VoxelBlock<FieldType> * block) {
          Eigen::Vector3i blockCoord = block->coordinates();
          unsigned int y, z, x; 
          Eigen::Vector3i blockSide = Eigen::Vector3i::Constant(se::VoxelBlock<FieldType>::side);
          Eigen::Vector3i start = blockCoord.cwiseMax(_min);
          Eigen::Vector3i last = (blockCoord + blockSide).cwiseMin(_max);

          for(z = start(2); z < last(2); ++z) {
            for (y = start(1); y < last(1); ++y) {
              for (x = start(0); x < last(0); ++x) {
                Eigen::Vector3i vox = Eigen::Vector3i(x, y, z);
                VoxelBlockHandler<FieldType> handler = {block, vox};
                _function(handler, vox);
              }
            }
          }
        }

        void update_node(se::Node<FieldType> * node) { 
          Eigen::Vector3i voxel = Eigen::Vector3i(unpack_morton(node->code));
#pragma omp simd
          for(int i = 0; i < 8; ++i) {
            const Eigen::Vector3i dir =  Eigen::Vector3i((i & 1) > 0, (i & 2) > 0, (i & 4) > 0);
            voxel = voxel + (dir * (node->side/2));
            if(!(se::math::in(voxel(0), _min(0), _max(0)) && 
                 se::math::in(voxel(1), _min(1), _max(1)) && 
                 se::math::in(voxel(2), _min(2), _max(2)))) continue;
            NodeHandler<FieldType> handler = {node, i};
            _function(handler, voxel);
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
        Eigen::Vector3i _min;
        Eigen::Vector3i _max;
      };

    /*!
     * \brief Applies a function object to each voxel/octant in the map. 
     * \param map Octree on which the function is going to be applied.
     * \param funct Update function to be applied.
     */
    template <typename FieldType, template <typename FieldT> class MapT, 
              typename UpdateF>
    void axis_aligned_map(MapT<FieldType>& map, UpdateF funct) {
    axis_aligned<FieldType, MapT, UpdateF> aa_functor(map, funct);
    aa_functor.apply();
    }

    template <typename FieldType, template <typename FieldT> class MapT, 
              typename UpdateF>
    void axis_aligned_map(MapT<FieldType>& map, UpdateF funct,
        const Eigen::Vector3i& min, const Eigen::Vector3i& max) {
    axis_aligned<FieldType, MapT, UpdateF> aa_functor(map, funct, min,  max);
    aa_functor.apply();
    }
  }
}
#endif
