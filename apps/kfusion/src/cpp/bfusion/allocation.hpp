#ifndef BFUSION_ALLOC_H
#define BFUSION_ALLOC_H
#include "math_utils.h"

/* Compute step size based on distance travelled along the ray */ 
static inline float compute_stepsize(const float dist_travelled, const float hf_band,
    const float voxelSize) {
  float new_step;
  float half_band = hf_band * 0.5f;
  if(dist_travelled < hf_band) new_step = voxelSize;
  else if(dist_travelled <= hf_band + half_band) new_step = 10.f * voxelSize; 
  else new_step = 30.f * voxelSize;
  return new_step;
}

/* Compute octree level given a step size */ 
static inline int step_to_depth(const float step, const int max_depth, 
    const float voxelsize) {
  return static_cast<int>(std::log2f(voxelsize/step)) + max_depth;
}
#endif
