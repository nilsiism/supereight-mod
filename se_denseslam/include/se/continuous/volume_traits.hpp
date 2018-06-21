#ifndef VOLUME_H
#define VOLUME_H

// Data types definitions
#include <se/voxel_traits.hpp>

/******************************************************************************
 *
 * KFusion Truncated Signed Distance Function voxel traits
 *
****************************************************************************/

class SDF {};
template<>
struct voxel_traits<SDF> {
  typedef float2 value_type;
  static inline value_type empty(){ return make_float2(1.f, -1.f); }
  static inline value_type initValue(){ return make_float2(1.f, 0.f); }
};

/******************************************************************************
 *
 * Bayesian Fusion voxel traits and algorithm specificic defines
 *
****************************************************************************/

class BFusion {};
template<>
struct voxel_traits<BFusion> {
  typedef struct  {
    float x;
    double y;
  } value_type;
  static inline value_type empty(){ return {0.f, 0.f}; }
  static inline value_type initValue(){ return {0.f, 0.f}; }
};

// Windowing parameters
#define DELTA_T   1.f
#define CAPITAL_T 4.f

#define INTERP_THRESH 0.05f
#define SURF_BOUNDARY 0.f
#define TOP_CLAMP     1000.f
#define BOTTOM_CLAMP  (-TOP_CLAMP)

#endif
