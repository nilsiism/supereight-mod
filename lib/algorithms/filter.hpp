#ifndef ACTIVE_LIST_HPP
#define ACTIVE_LIST_HPP

#include <math_utils.h> 
#include <node.hpp>
#include <memory_pool.hpp>
#include <utils/morton_utils.hpp>

namespace algorithms {

  template <typename VoxelBlockType>
    bool in_frustum(const VoxelBlockType* v, float voxelSize, 
        const Matrix4& camera, const int2& frameSize) {
      const float3 block_coord = make_float3(v->coordinates()) * voxelSize;
      const float3 v_camera = camera*block_coord;
      const int2 px = make_int2(v_camera.x/v_camera.z, v_camera.y/v_camera.z);
      if( px.x >= 0 && px.x < frameSize.x && px.y >= 0 && px.y < frameSize.y)
        return true;
      return false;
    }

  template <typename ValueType, typename P>
    bool satisfies(const ValueType& el, P predicate) {
      return predicate(el);
    }

  template <typename ValueType, typename P, typename... Ps>
    bool satisfies(const ValueType& el, P predicate, Ps... others) {
      return predicate(el) || satisfies(el, others...);
    }

  template <typename BlockType, typename... Predicates>
    void filter(std::vector<BlockType *>& out,
        const MemoryPool<BlockType>& block_array, Predicates... ps) {
      for(unsigned int i = 0; i < block_array.size(); ++i) {
        if(satisfies(block_array[i], ps...)){
          out.push_back(block_array[i]);
        }
      } 
    }
}
#endif
