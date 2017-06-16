#ifndef VOLUME_H
#define VOLUME_H

// Continuous volume template implementation
#include "volume_template.hpp"

// Discrete octree implementation
#include "octree.hpp"

// Data types definitions
#include "voxel_traits.hpp"

class SDF;
template<>
struct voxel_traits<SDF> {
  typedef float2 ComputeType;
  typedef short2 StoredType;
  static inline ComputeType empty(){ return make_float2(1.f, -1.f); }
  static inline StoredType initValue(){ return make_short2(32766, 0); }
  static inline StoredType translate(const ComputeType value) {
     return make_short2(value.x * 32766.0f, value.y);
  }
  static inline ComputeType translate(const StoredType value){
    return make_float2(value.x * 0.00003051944088f, value.y);
  }
};

class BFusion;
template<>
struct voxel_traits<BFusion> {
  typedef float1 ComputeType;
  typedef float1 StoredType;
  static inline ComputeType empty(){ return make_float1(0.f); }
  static inline StoredType initValue(){ return make_float1(0.f); }
  static inline StoredType translate(const ComputeType value) {
     return value;
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
