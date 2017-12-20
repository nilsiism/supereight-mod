#ifndef VOLUME_TEMPLATE_H
#define VOLUME_TEMPLATE_H

#include <voxel_traits.hpp>
#include <memory_pool.hpp>
#include <data.hpp>
#include <type_traits>
#include <algorithms/alloc_list.hpp>
#include <algorithms/filter.hpp>
#include <algorithms/mapping.hpp>
#include <functors/projective_functor.hpp>
#include <functional>
#include <cstring>
#include "../mapping.hpp"
#include "../bfusion/allocation.hpp"

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

    compute_type get(const float3 & p) const {
      const float inverseVoxelSize = _size/_dim;
      const int3 scaled_pos = make_int3(make_float3((p.x * inverseVoxelSize),
          (p.y * inverseVoxelSize), (p.z * inverseVoxelSize)));
      return _map_index.get_fine(scaled_pos.x, scaled_pos.y, scaled_pos.z);
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
                            const int ) {

      using namespace std::placeholders;
      int num_vox_per_pix = _dim/((VoxelBlock<FieldType>::side)*(_dim/_size));
      

      size_t total = num_vox_per_pix * frameSize.x * frameSize.y;
      _allocationList.reserve(total);
      static bool initialised = false;
      if(!initialised) {
        std::memset(_allocationList.data(), 0, sizeof(unsigned int) * total);
        initialised = true;
      }

      const float hf_band = 2*mu;

      int allocated = buildOctantList(_allocationList.data(),
          _allocationList.capacity(), _map_index, pose, K, depthmap, frameSize,
          _dim/_size, compute_stepsize, step_to_depth,  hf_band);
      _map_index.alloc_update(_allocationList.data(), allocated);
      
      // double current = frame * (1.f/30.f);
      // struct bfusion_update funct(depthmap, frameSize, mu, current);
      // iterators::projective_functor<FieldType, Indexer, struct bfusion_update> 
      //   it(_map_index, funct, inverse(pose), K, make_int2(frameSize));
      // it.apply();

      struct sdf_update funct(depthmap, frameSize, mu, 100);
      iterators::projective_functor<FieldType, Indexer, struct sdf_update> 
        it(_map_index, funct, inverse(pose), K, make_int2(frameSize));
      it.apply();
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
