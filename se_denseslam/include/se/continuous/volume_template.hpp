#ifndef VOLUME_TEMPLATE_H
#define VOLUME_TEMPLATE_H

#include <se/voxel_traits.hpp>
#include <se/utils/memory_pool.hpp>
#include <type_traits>
#include <cstring>

template <typename T>
class Void {};

/**
 * Continuous volume abstraction
 * Sparse, dynamically allocated storage accessed through the 
 * appropriate indexer (octree/hash table).
 * */ 
template <typename FieldType, template<typename> class Indexer> 
class VolumeTemplate {

  public:
    typedef voxel_traits<FieldType> traits_type;
    typedef typename traits_type::value_type value_type;
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

    void set(const uint3 & , const value_type& ) {}

    value_type operator[](const float3 & p) const {
      const float inverseVoxelSize = _size/_dim;
      const int3 scaled_pos = make_int3(make_float3((p.x * inverseVoxelSize),
          (p.y * inverseVoxelSize), (p.z * inverseVoxelSize)));
      return _map_index.get(scaled_pos.x, scaled_pos.y, scaled_pos.z);
    }

    value_type get(const float3 & p) const {
      const float inverseVoxelSize = _size/_dim;
      const int3 scaled_pos = make_int3(make_float3((p.x * inverseVoxelSize),
          (p.y * inverseVoxelSize), (p.z * inverseVoxelSize)));
      return _map_index.get_fine(scaled_pos.x, scaled_pos.y, scaled_pos.z);
    }

    value_type operator[](const uint3 p) const {
      return _map_index.get(p.x, p.y, p.z);
    }

    template <typename FieldSelector>
    float interp(const float3 & pos, FieldSelector select) const {
      const float inverseVoxelSize = _size / _dim;
      const Eigen::Vector3f scaled_pos((pos.x * inverseVoxelSize),
          (pos.y * inverseVoxelSize),
          (pos.z * inverseVoxelSize));
      return _map_index.interp(scaled_pos, select);
    }

    template <typename FieldSelector>
    Eigen::Vector3f grad(const float3 & pos, FieldSelector select) const {

      const float inverseVoxelSize = _size / _dim;
      const Eigen::Vector3f scaled_pos((pos.x * inverseVoxelSize),
          (pos.y * inverseVoxelSize),
          (pos.z * inverseVoxelSize));
      return _map_index.grad(scaled_pos, select);
    }

    unsigned int _size;
    float _dim;
    std::vector<se::key_t> _allocationList;
    Indexer<FieldType> _map_index; 

  private:

    inline uint3 pos(const float3 & p) const {
      static const float inverseVoxelSize = _size/_dim;
      return make_uint3(p.x * inverseVoxelSize, p.y * inverseVoxelSize, 
          p.z * inverseVoxelSize);
    }
};
#endif