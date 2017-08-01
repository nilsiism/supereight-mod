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

#ifdef _OPENMP
  template <typename BlockType, typename... Predicates>
    void filter(std::vector<BlockType *>& out,
        const MemoryPool<BlockType>& block_array, Predicates... ps) {

      std::vector<BlockType *> temp;
      int num_elem = block_array.size();
      temp.resize(num_elem);

      int * thread_start = new int[omp_get_max_threads()];
      int * thread_end = new int[omp_get_max_threads()];
      int spawn_threads;
#pragma omp parallel
      {
        int threadid = omp_get_thread_num(); 
        int num_threads = omp_get_num_threads();
        int my_start = thread_start[threadid] = (threadid) * num_elem / num_threads;
        int my_end   = (threadid+1) * num_elem / num_threads;
        int count = 0;
        for(int i = my_start; i < my_end; ++i) {
          if(satisfies(block_array[i], ps...)){
            temp[my_start + count] = block_array[i];
            count++;
          }
        } 
        /* Store the actual end */
        thread_end[threadid] = count;
        if(threadid == 0) spawn_threads = num_threads;
      }
      
      int total = 0;
      for(int i = 0; i < spawn_threads; ++i) {
        total += thread_end[i];
      }
      out.resize(total);
      /* Copy the first */
      std::memcpy(out.data(), temp.data(), sizeof(BlockType *) * thread_end[0]);
      int copied = thread_end[0];
      /* Copy the rest */
      for(int i = 1; i < spawn_threads; ++i) {
        std::memcpy(out.data() + copied, 
            temp.data() + thread_start[i], sizeof(BlockType *) * thread_end[i]);
        copied += thread_end[i];
      }
    }

#else
  template <typename BlockType, typename... Predicates>
    void filter(std::vector<BlockType *>& out,
        const MemoryPool<BlockType>& block_array, Predicates... ps) {
      for(unsigned int i = 0; i < block_array.size(); ++i) {
        if(satisfies(block_array[i], ps...)){
          out.push_back(block_array[i]);
        }
      } 
    }
#endif
}
#endif
