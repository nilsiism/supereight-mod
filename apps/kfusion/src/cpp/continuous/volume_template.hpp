#ifndef VOLUME_TEMPLATE_H
#define VOLUME_TEMPLATE_H

#include <voxel_traits.hpp>
#include <memory_pool.hpp>
#include <data.hpp>
#include <type_traits>
#include <algorithms/alloc_list.hpp>
#include <algorithms/filter.hpp>
// #include <algorithms/mapping.hpp>
#include <functional>
#include "../mapping.hpp"

template <typename T>
class Void {};

template <typename FieldType, typename StorageLayout, template<typename> class Indexer>
class VolumeTemplate {};

/**
 * Continuous volume abstraction
 * Sparse, dynamically allocated storage accessed through the 
 * appropriate indexer (octree/hash table).
 * */ 
template <typename FieldType, template<typename> class Indexer> 
class VolumeTemplate<FieldType, DynamicStorage, Indexer> {

  public:
    typedef voxel_traits<FieldType> traits_type;
    typedef typename traits_type::ComputeType compute_type;
    typedef typename traits_type::StoredType stored_type;
    typedef FieldType field_type;

    void init(uint s, float d) {
      _size = s;
      _dim = d; 
      _map_index.init(_size, _dim);
    }

    void release(){};

    inline float3 pos(const uint3 & p) const {
      static const float voxelSize = _dim/_size;
      return make_float3(p.x * voxelSize, p.y * voxelSize, p.z * voxelSize);
    }

    void set(const uint3 & pos, const compute_type& d) {}

    compute_type operator[](const float3 & p) const {
      const float inverseVoxelSize = _size/_dim;
      const int3 scaled_pos = make_int3(make_float3((p.x * inverseVoxelSize),
          (p.y * inverseVoxelSize), (p.z * inverseVoxelSize)));
      return _map_index.get(scaled_pos.x, scaled_pos.y, scaled_pos.z);
    }

    compute_type operator[](const uint3 p) const {
      return _map_index.get(p.x, p.y, p.z);
    }

    template <typename FieldSelector>
    float interp(const float3 & pos, FieldSelector select) const {
      const float inverseVoxelSize = _size / _dim;
      const float3 scaled_pos = make_float3((pos.x * inverseVoxelSize),
          (pos.y * inverseVoxelSize),
          (pos.z * inverseVoxelSize));
      return _map_index.interp(scaled_pos, select);
    }

    template <typename FieldSelector>
    float3 grad(const float3 & pos, FieldSelector select) const {

      const float inverseVoxelSize = _size / _dim;
      const float3 scaled_pos = make_float3((pos.x * inverseVoxelSize),
          (pos.y * inverseVoxelSize),
          (pos.z * inverseVoxelSize));
      return _map_index.grad(scaled_pos, select);
    }

    void updateVolume(const Matrix4 &pose, const Matrix4& K,
                            const float *depthmap,
                            const uint2 frameSize, const float mu,
                            const int frame) {

      using namespace std::placeholders;
      int num_vox_per_pix = (3 * mu)/(_dim/_size);

     auto compute_sdf = std::bind(integrate_bfusion, _1,  depthmap, frameSize,
         _dim/_size, inverse(pose), K,  mu, 100);

      _allocationList.reserve(num_vox_per_pix * frameSize.x * frameSize.y);
      const int allocated = 
        buildAllocationList(_allocationList.data(), _allocationList.capacity(),  
          _map_index, pose, K, depthmap, frameSize, _size, 
          _dim/_size, mu);
      _map_index.alloc_update(_allocationList.data(), allocated, 7, compute_sdf);

     const int max_depth = 5;
     const float step = _dim/std::pow(2, 5);
     const int updated = 
       buildIntersectionList(_allocationList.data(), _allocationList.capacity(),  
         _map_index, pose, K, depthmap, frameSize, max_depth, 
         _dim/_size, step, mu);
     _map_index.alloc_update(_allocationList.data(), updated, max_depth, compute_sdf);

      std::vector<VoxelBlock<FieldType> *> active_list;
      const MemoryPool<VoxelBlock<FieldType> >& block_array = 
        _map_index.getBlockBuffer();

      auto in_frustum_predicate = 
        std::bind(algorithms::in_frustum<VoxelBlock<FieldType>>, _1, 
         _dim/_size, K*inverse(pose), make_int2(frameSize)); 
      auto is_active_predicate = [](const VoxelBlock<FieldType>* b) {
        return b->active();
      };
      
      algorithms::filter(active_list, block_array, is_active_predicate, in_frustum_predicate);
      VoxelBlock<FieldType> ** list = active_list.data();
      unsigned int num_active = active_list.size();
      integratePass(list, num_active, depthmap, frameSize, _dim/_size,
          inverse(pose), K,  mu, 100, frame);

      MemoryPool<Node<FieldType> >& nodes_array =
        _map_index.getNodesBuffer();
      const int size = nodes_array.size();
      integratePass(nodes_array, size, compute_sdf);
    }

    unsigned int _size;
    float _dim;
    std::vector<uint> _allocationList;
    Indexer<FieldType> _map_index; 

  private:

    inline uint3 pos(const float3 & p) const {
      static const float inverseVoxelSize = _size/_dim;
      return make_uint3(p.x * inverseVoxelSize, p.y * inverseVoxelSize, 
          p.z * inverseVoxelSize);
    }
};
#endif
