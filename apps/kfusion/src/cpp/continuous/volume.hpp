#ifndef VOLUME_H
#define VOLUME_H

// Continuous volume template implementation
#include "volume_template.hpp"

// Discrete octree implementation
#include <octree.hpp>

// Data types definitions
#include "voxel_traits.hpp"

/******************************************************************************
 *
 * KFusion Truncated Signed Distance Function voxel traits
 *
****************************************************************************/

class SDF {};
template<>
struct voxel_traits<SDF> {
  typedef float2 ComputeType;
  typedef float2 StoredType;
  static inline ComputeType empty(){ return make_float2(1.f, -1.f); }
  static inline StoredType initValue(){ return make_float2(1.f, 0.f); }
  static inline ComputeType translate(const StoredType value){
    return make_float2(value.x, value.y);
  }
};

/******************************************************************************
 *
 * Bayesian Fusion voxel traits and algorithm specificic defines
 *
****************************************************************************/

class BFusion {};
template<>
struct voxel_traits<BFusion> {
  typedef float2 ComputeType;
  typedef float2 StoredType;
  static inline ComputeType empty(){ return make_float2(-101.f, 0.f); }
  static inline StoredType initValue(){ return make_float2(0.f, 0.f); }
  static inline StoredType translate(const ComputeType value) {
     return value;
  }
};

// Windowing parameters
#define DELTA_T   1.f
#define CAPITAL_T 50.f

#define INTERP_THRESH 0.05f
#define SURF_BOUNDARY 0.f
#define TOP_CLAMP     100.f
#define BOTTOM_CLAMP  (-TOP_CLAMP)

/******************************************************************************
 *
 * Default typedefs for the Volume Class
 * 
******************************************************************************/


#ifndef STORAGE_LAYOUT
#define STORAGE_LAYOUT DynamicStorage
#endif

#ifndef INDEX_STRUCTURE
#define INDEX_STRUCTURE Octree
#endif

typedef BFusion FieldType;

template <typename T>
using Volume = VolumeTemplate<T, STORAGE_LAYOUT, INDEX_STRUCTURE>;
#endif
