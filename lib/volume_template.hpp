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
 * Contiguous, preallocated storage layout specified by the volume
 * discretisation resolution
 */ 
template <typename FieldType > 
class VolumeTemplate<FieldType, StaticStorage, Void> {

  public:
    typedef kfusion_voxel_traits<FieldType> traits_type;
    typedef typename traits_type::ComputeType compute_type;
    typedef typename traits_type::StoredType stored_type;

    void init(uint s, float d) {
      _size = s;
      _dim = d;
      _data.init(s * s * s);
      for (unsigned int x = 0; x < _size; x++)
        for (unsigned int y = 0; y < _size; y++) {
          for (unsigned int z = 0; z < _size; z++) {
            stored_type& d = _data(x + y*_size + z*_size*_size); 
            d = traits_type::initValue(); 
          }
        }
    }

    void release(){ 
      _data.release();
    };

    inline float3 pos(const uint3 & p) const {
      static const float voxelSize = _dim/_size;
      return make_float3((p.x + 0.5f) * voxelSize, (p.y + 0.5f) * voxelSize, 
          (p.z + 0.5f) * voxelSize);
    }

    void set(const uint3 & pos, const compute_type& d) {
      stored_type& data = _data(pos.x + pos.y * _size + pos.z * _size * _size);
      data = traits_type::translate(d);
    }

    compute_type operator[](const float3 & p) const {
      const float3 scaled_pos = make_float3((p.x * _size / _dim),
          (p.y * _size / _dim),
          (p.z * _size / _dim));
      const int3 pos = max(make_int3(floorf(scaled_pos)), make_int3(0));
      const stored_type d = _data(pos.x + pos.y * _size + pos.z * _size * _size);
      return traits_type::translate(d);
    }

    compute_type operator[](const uint3 & pos) const {
      const stored_type d = _data(pos.x + pos.y * _size + pos.z * _size * _size);
      return traits_type::translate(d); 
    }

    float interp(const float3 & pos) const {
      const float inverseVoxelSize = _size / _dim;
      const float3 scaled_pos = make_float3((pos.x * inverseVoxelSize) - 0.5f,
          (pos.y * inverseVoxelSize) - 0.5f,
          (pos.z * inverseVoxelSize) - 0.5f);
      const int3 base = make_int3(floorf(scaled_pos));
      const float3 factor = fracf(scaled_pos);
      const int3 lower = max(base, make_int3(0));
      const int3 upper = min(base + make_int3(1),
          make_int3(_size) - make_int3(1));
      return (((data_raw(lower.x, lower.y, lower.z).x * (1 - factor.x)
              + data_raw(upper.x, lower.y, lower.z).x * factor.x) * (1 - factor.y)
            + (data_raw(lower.x, upper.y, lower.z).x * (1 - factor.x)
              + data_raw(upper.x, upper.y, lower.z).x * factor.x) * factor.y)
          * (1 - factor.z)
          + ((data_raw(lower.x, lower.y, upper.z).x * (1 - factor.x)
              + data_raw(upper.x, lower.y, upper.z).x * factor.x)
            * (1 - factor.y)
            + (data_raw(lower.x, upper.y, upper.z).x * (1 - factor.x)
              + data_raw(upper.x, upper.y, upper.z).x * factor.x)
            * factor.y) * factor.z) * 0.00003051944088f;

    }

    float3 grad(const float3 & pos) const {

      const float voxelSize = _dim / _size;
      const float inverseVoxelSize = _size / _dim;
      const float3 scaled_pos = make_float3((pos.x * inverseVoxelSize) - 0.5f,
          (pos.y * inverseVoxelSize) - 0.5f,
          (pos.z * inverseVoxelSize) - 0.5f);
      const int3 base = make_int3(floorf(scaled_pos));
      const float3 factor = fracf(scaled_pos);
      const int3 lower_lower = max(base - make_int3(1), make_int3(0));
      const int3 lower_upper = max(base, make_int3(0));
      const int3 upper_lower = min(base + make_int3(1),
          make_int3(_size) - make_int3(1));
      const int3 upper_upper = min(base + make_int3(2),
          make_int3(_size) - make_int3(1));
      const int3 & lower = lower_upper;
      const int3 & upper = upper_lower;

      float3 gradient;

      gradient.x = (((data_raw(upper_lower.x, lower.y, lower.z).x
              - data_raw(lower_lower.x, lower.y, lower.z).x) * (1 - factor.x)
            + (data_raw(upper_upper.x, lower.y, lower.z).x
              - data_raw(lower_upper.x, lower.y, lower.z).x) * factor.x)
          * (1 - factor.y)
          + ((data_raw(upper_lower.x, upper.y, lower.z).x
              - data_raw(lower_lower.x, upper.y, lower.z).x) * (1 - factor.x)
            + (data_raw(upper_upper.x, upper.y, lower.z).x
              - data_raw(lower_upper.x, upper.y, lower.z).x)
            * factor.x) * factor.y) * (1 - factor.z)
        + (((data_raw(upper_lower.x, lower.y, upper.z).x
                - data_raw(lower_lower.x, lower.y, upper.z).x) * (1 - factor.x)
              + (data_raw(upper_upper.x, lower.y, upper.z).x
                - data_raw(lower_upper.x, lower.y, upper.z).x)
              * factor.x) * (1 - factor.y)
            + ((data_raw(upper_lower.x, upper.y, upper.z).x
                - data_raw(lower_lower.x, upper.y, upper.z).x)
              * (1 - factor.x)
              + (data_raw(upper_upper.x, upper.y, upper.z).x
                - data_raw(lower_upper.x, upper.y, upper.z).x)
              * factor.x) * factor.y) * factor.z;

      gradient.y = (((data_raw(lower.x, upper_lower.y, lower.z).x
              - data_raw(lower.x, lower_lower.y, lower.z).x) * (1 - factor.x)
            + (data_raw(upper.x, upper_lower.y, lower.z).x
              - data_raw(upper.x, lower_lower.y, lower.z).x) * factor.x)
          * (1 - factor.y)
          + ((data_raw(lower.x, upper_upper.y, lower.z).x
              - data_raw(lower.x, lower_upper.y, lower.z).x) * (1 - factor.x)
            + (data_raw(upper.x, upper_upper.y, lower.z).x
              - data_raw(upper.x, lower_upper.y, lower.z).x)
            * factor.x) * factor.y) * (1 - factor.z)
        + (((data_raw(lower.x, upper_lower.y, upper.z).x
                - data_raw(lower.x, lower_lower.y, upper.z).x) * (1 - factor.x)
              + (data_raw(upper.x, upper_lower.y, upper.z).x
                - data_raw(upper.x, lower_lower.y, upper.z).x)
              * factor.x) * (1 - factor.y)
            + ((data_raw(lower.x, upper_upper.y, upper.z).x
                - data_raw(lower.x, lower_upper.y, upper.z).x)
              * (1 - factor.x)
              + (data_raw(upper.x, upper_upper.y, upper.z).x
                - data_raw(upper.x, lower_upper.y, upper.z).x)
              * factor.x) * factor.y) * factor.z;

      gradient.z = (((data_raw(lower.x, lower.y, upper_lower.z).x
              - data_raw(lower.x, lower.y, lower_lower.z).x) * (1 - factor.x)
            + (data_raw(upper.x, lower.y, upper_lower.z).x
              - data_raw(upper.x, lower.y, lower_lower.z).x) * factor.x)
          * (1 - factor.y)
          + ((data_raw(lower.x, upper.y, upper_lower.z).x
              - data_raw(lower.x, upper.y, lower_lower.z).x) * (1 - factor.x)
            + (data_raw(upper.x, upper.y, upper_lower.z).x
              - data_raw(upper.x, upper.y, lower_lower.z).x)
            * factor.x) * factor.y) * (1 - factor.z)
        + (((data_raw(lower.x, lower.y, upper_upper.z).x
                - data_raw(lower.x, lower.y, lower_upper.z).x) * (1 - factor.x)
              + (data_raw(upper.x, lower.y, upper_upper.z).x
                - data_raw(upper.x, lower.y, lower_upper.z).x)
              * factor.x) * (1 - factor.y)
            + ((data_raw(lower.x, upper.y, upper_upper.z).x
                - data_raw(lower.x, upper.y, lower_upper.z).x)
              * (1 - factor.x)
              + (data_raw(upper.x, upper.y, upper_upper.z).x
                - data_raw(upper.x, upper.y, lower_upper.z).x)
              * factor.x) * factor.y) * factor.z;

      return gradient
        * make_float3(voxelSize)
        * (0.5f * 0.00003051944088f);
    }

    unsigned int _size;
    float _dim;
    Data<FieldType, StaticStorage> _data; // This should actually be managed by the Memory Pool.

  private:

    inline uint3 pos(const float3 & p) const {
      static const float inverseVoxelSize = _size/_dim;
      return make_uint3(p.x * inverseVoxelSize, p.y * inverseVoxelSize, 
          p.z * inverseVoxelSize);
    }

    stored_type data_raw(const uint x, const uint y, const uint z) const {
      const stored_type d = _data(x + y * _size + z * _size * _size);
      return d; 
    }

};

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
      _map_index.init(_size, make_float3(_dim));
    }

    void release(){};

    inline float3 pos(const uint3 & p) const {
      static const float voxelSize = _dim/_size;
      return make_float3((p.x + 0.5f) * voxelSize, (p.y + 0.5f) * voxelSize, 
          (p.z + 0.5f) * voxelSize);
    }

    void set(const uint3 & pos, const compute_type& d) {}

    compute_type operator[](const float3 & p) const {
      const float3 scaled_pos = make_float3((p.x * _size / _dim),
          (p.y * _size / _dim),
          (p.z * _size / _dim));
      return _map_index.get(scaled_pos);
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

    void updateVolume(const Matrix4 &pose, const float4& k,
                            const float *depthmap,
                            const uint2 frameSize, const float mu,
                            const int frame) {

      int num_vox_per_pix = (2 * mu)/(_dim/_size);
      _allocationList.reserve(num_vox_per_pix * frameSize.x * frameSize.y);
      const int allocated = 
        buildAllocationList(_allocationList.data(), _allocationList.capacity(),  
          _map_index, pose, getCameraMatrix(k), depthmap, frameSize, _size, 
          _dim/_size, mu, frame);
      _map_index.allocate(_allocationList.data(), allocated);
      _map_index.integrateFrame(pose, k, depthmap, frameSize, mu, frame); 
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
