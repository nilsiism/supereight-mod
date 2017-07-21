#ifndef VOLUME_H
#define VOLUME_H

// Continuous volume template implementation
#include "volume_template.hpp"

// Discrete octree implementation
#include "octree.hpp"

// Data types definitions
#include "voxel_traits.hpp"

template<>
struct voxel_traits<SDF> {
  typedef float2 ComputeType;
  typedef float2 StoredType;
  static inline ComputeType empty(){ return make_float2(1.f, -1.f); }
  static inline StoredType initValue(){ return make_float2(1.f, 0); }
  static inline StoredType translate(const ComputeType value) {
     return make_float2(value.x, value.y);
  }
};

/**
 * Default typedefs for the Volume Class
 * */

#ifndef FIELD_TYPE
#define FIELD_TYPE SDF
#endif

#ifndef STORAGE_LAYOUT
#define STORAGE_LAYOUT DynamicStorage
#endif

#ifndef INDEX_STRUCTURE
#define INDEX_STRUCTURE Octree
#endif

typedef VolumeTemplate<FIELD_TYPE, STORAGE_LAYOUT, INDEX_STRUCTURE> Volume;
#endif