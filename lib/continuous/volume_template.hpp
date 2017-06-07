#ifndef VOLUME_TEMPLATE_H
#define VOLUME_TEMPLATE_H

#include <voxel_traits.hpp>
#include <memory_pool.hpp>
#include <data.hpp>
#include <type_traits>
#include <algorithms/alloc_list.hpp>

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
    typedef kfusion_voxel_traits<FieldType> traits_type;
    typedef typename traits_type::ComputeType compute_type;
    typedef typename traits_type::StoredType stored_type;

    void init(uint s, float d) {
      _size = s;
      _dim = d; 
      _map_index.init(_size, _dim);
    }

    void release(){};

    inline float3 pos(const uint3 & p) const {
      static const float voxelSize = _dim/_size;
      return make_float3((p.x + 0.5f) * voxelSize, (p.y + 0.5f) * voxelSize, 
          (p.z + 0.5f) * voxelSize);
    }

    void set(const uint3 & pos, const compute_type& d) {}

    compute_type operator[](const float3 & p) const {
      const int3 scaled_pos = make_int3(make_float3((p.x * _size / _dim),
          (p.y * _size / _dim),
          (p.z * _size / _dim)));
      return _map_index.get(scaled_pos.x, scaled_pos.y, scaled_pos.z);
    }

    compute_type operator[](const uint3 p) const {
      return _map_index.get(p.x, p.y, p.z);
    }

    float interp(const float3 & pos) const {
      const float inverseVoxelSize = _size / _dim;
      const float3 scaled_pos = make_float3((pos.x * inverseVoxelSize) - 0.5f,
          (pos.y * inverseVoxelSize) - 0.5f,
          (pos.z * inverseVoxelSize) - 0.5f);
      return _map_index.interp(scaled_pos);
    }

    float3 grad(const float3 & pos) const {

      const float voxelSize = _dim / _size;
      const float inverseVoxelSize = _size / _dim;
      const float3 scaled_pos = make_float3((pos.x * inverseVoxelSize) - 0.5f,
          (pos.y * inverseVoxelSize) - 0.5f,
          (pos.z * inverseVoxelSize) - 0.5f);
      return _map_index.grad(scaled_pos);
    }

    void updateVolume(const Matrix4 &pose, const Matrix4& K,
                            const float *depthmap,
                            const uint2 frameSize, const float mu,
                            const int frame) {

      int num_vox_per_pix = (2 * mu)/(_dim/_size);
      _allocationList.reserve(num_vox_per_pix * frameSize.x * frameSize.y);
      const int allocated = 
        buildAllocationList(_allocationList.data(), _allocationList.capacity(),  
          _map_index, pose, K, depthmap, frameSize, _size, 
          _dim/_size, mu, frame);
      _map_index.allocate(_allocationList.data(), allocated);
      _map_index.integrateFrame(pose, K, depthmap, frameSize, mu, frame); 
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
