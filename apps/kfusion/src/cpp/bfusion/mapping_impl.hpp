#ifndef BFUSION_MAPPING_HPP
#define BFUSION_MAPPING_HPP

#include <node.hpp>
#include <constant_parameters.h>
#include "bspline_lookup.cc"
#include "../continuous/volume_traits.hpp"

float interpDepth(const float * depth, const uint2 depthSize, 
    const float2 proj) {
  // https://en.wikipedia.org/wiki/Bilinear_interpolation

  // Pixels version
  const float x1 = (floorf(proj.x));
  const float y1 = (floorf(proj.y + 1));
  const float x2 = (floorf(proj.x + 1));
  const float y2 = (floorf(proj.y));

  // Half pixels
  // const float x1 = (float) (int(proj.x - 0.5f)) + 0.5f;
  // const float y1 = (float) (int(proj.y + 0.5f)) + 0.5f;
  // const float x2 = (float) (int(proj.x + 0.5f)) + 0.5f;
  // const float y2 = (float) (int(proj.y - 0.5f)) + 0.5f;

  const float d11 = depth[int(x1) +  depthSize.x*int(y1)];
  const float d12 = depth[int(x1) +  depthSize.x*int(y2)];
  const float d21 = depth[int(x2) +  depthSize.x*int(y1)];
  const float d22 = depth[int(x2) +  depthSize.x*int(y2)];
  
  if( d11 == 0.f || d12 == 0.f || d21 == 0.f || d22 == 0.f ) return 0.f;

  const float f11 = 1.f / d11;
  const float f12 = 1.f / d12;
  const float f21 = 1.f / d21;
  const float f22 = 1.f / d22;
  
  // Filtering version
  const float d =  1.f / 
                    ( (   f11 * (x2 - proj.x) * (y2 - proj.y)
                        + f21 * (proj.x - x1) * (y2 - proj.y)
                        + f12 * (x2 - proj.x) * (proj.y - y1)
                        + f22 * (proj.x - x1) * (proj.y - y1)
                      ) / ((x2 - x1) * (y2 - y1))
                    );

  static const float interp_thresh = 0.05f;
  if (fabs(d - d11) < interp_thresh && fabs(d - d12) < interp_thresh &&
      fabs(d - d21) < interp_thresh && fabs(d - d22) < interp_thresh) 
    return d;
  else 
    return depth[int(proj.x + 0.5f) + depthSize.x*int(proj.y+0.5f)];

  // Non-Filtering version
  // return  1.f / 
  //         ( (   f11 * (x2 - proj.x) * (y2 - proj.y)
  //             + f21 * (proj.x - x1) * (y2 - proj.y)
  //             + f12 * (x2 - proj.x) * (proj.y - y1)
  //             + f22 * (proj.x - x1) * (proj.y - y1)
  //           ) / ((x2 - x1) * (y2 - y1))
  //         );
}

static inline float bspline(float t){
  float value = 0.f;
  if(t >= -3.0f && t <= -1.0f) {
    value = std::pow((3 + t), 3)/48.0f;   
  } else if( t > -1 && t <= 1) {
    value = 0.5f + (t*(3 + t)*(3 - t))/24.f;
  } else if(t > 1 && t <= 3){
    value = 1 - std::pow((3 - t), 3)/48.f;
  } else if(t > 3) {
    value = 1.f;
  }
  return value;
}

static inline float H(const float val){
  const float Q_1 = bspline(val);
  const float Q_2 = bspline(val - 3);
  return Q_1 - Q_2 * 0.5f;
}

static const double const_offset =  0.0000001f;
const float scale_factor = (1.f - (farPlane - nearPlane) * const_offset); 

static inline float const_offset_integral(float t){
  float value = 0.f;
  if (nearPlane <= t && t <= farPlane)
      return (t - nearPlane) * const_offset;
  else if (farPlane < t)
      return (farPlane - nearPlane) * const_offset;
  return value;
}

static inline float __device__ bspline_memoized(float t){
  float value = 0.f;
  constexpr float inverseRange = 1/6.f;
  if(t >= -3.0f && t <= 3.0f) {
    unsigned int idx = ((t + 3.f)*inverseRange)*(bspline_num_samples - 1) + 0.5f;
    return bspline_lookup[idx];
  } 
  else if(t > 3) {
    value = 1.f;
  }
  return value;
}

static inline float HNew(const float val,const  float d_xr){
  const float Q_1 = bspline_memoized(val)    * scale_factor + const_offset_integral(d_xr      );
  const float Q_2 = bspline_memoized(val - 3)* scale_factor + const_offset_integral(d_xr - 3.f);
  return Q_1 - Q_2 * 0.5f;
}

static inline float updateLogs(const float prior, const float sample){
  return (prior + log2(sample / (1.f - sample)));
}

static inline float applyWindow(const float occupancy, const float boundary, 
    const float delta_t, const float tau){
  float fraction = 1.f / (1.f + delta_t / tau);
  return occupancy * fraction + boundary * (1.f - fraction);
}

void integrate(VoxelBlock<BFusion> * block, const float * depth, uint2 depthSize, 
    const float voxelSize, const Matrix4 world2cam, const Matrix4 K, 
    const float noiseFactor, const float timestamp) {

  const float3 delta = rotate(world2cam, make_float3(voxelSize, 0, 0));
  const float3 cameraDelta = rotate(K, delta);
  const int3 blockCoord = block->coordinates();
  bool is_visible = false;
  block->active(is_visible);

  unsigned int y, z, blockSide; 
  blockSide = VoxelBlock<BFusion>::side;
  unsigned int ylast = blockCoord.y + blockSide;
  unsigned int zlast = blockCoord.z + blockSide;
  for(z = blockCoord.z; z < zlast; ++z)
    for (y = blockCoord.y; y < ylast; ++y){
      int3 pix = make_int3(blockCoord.x, y, z);
      float3 start = world2cam * make_float3((pix.x) * voxelSize, 
          (pix.y) * voxelSize, (pix.z) * voxelSize);
      float3 camerastart = K * start;
#pragma omp simd
      for (unsigned int x = 0; x < blockSide; ++x){
        pix.x = x + blockCoord.x; 
        const float3 camera_voxel = camerastart + (x*cameraDelta);
        const float3 pos = start + (x*delta);
        if (pos.z < 0.0001f) continue;

        // Interpolation version
        const float2 pixel = make_float2( camera_voxel.x / camera_voxel.z + 0.5f,
            camera_voxel.y / camera_voxel.z + 0.5f);
        if (pixel.x < 0.5f || pixel.x > depthSize.x - 1.5f || 
            pixel.y < 0.5f || pixel.y > depthSize.y - 1.5f) continue;
        // const float depthSample = interpDepth(depth, depthSize, pixel);
        const uint2 px = make_uint2(pixel.x, pixel.y);
        const float depthSample = depth[px.x + depthSize.x*px.y];
        if (depthSample <=  0) continue;

        const float diff = (camera_voxel.z - depthSample)
          * std::sqrt( 1 + sq(pos.x / pos.z) + sq(pos.y / pos.z));
        const float sample = HNew(diff/(noiseFactor*sq(camera_voxel.z)), camera_voxel.z);
        if(sample == 0.5f) continue;
        is_visible = true;
        typename VoxelBlock<BFusion>::compute_type data = block->data(pix); // should do an offsetted access here
        data.x = applyWindow(data.x, SURF_BOUNDARY, DELTA_T, CAPITAL_T);
        data.x = clamp(updateLogs(data.x, sample), BOTTOM_CLAMP, TOP_CLAMP);
        data.y = timestamp;
        block->data(pix, data);
      }
    }
  block->active(is_visible);
}

void integrate_bfusion(Node<BFusion> * node, const int x, const int y, const int z,
    const float * depth, uint2 depthSize, 
    const float voxelSize, const Matrix4 invTrack, const Matrix4 K, 
    const float noiseFactor, const float ) { 

  const int3 voxel = make_int3(x, y, z);
  float3 pos = invTrack * (make_float3(voxel) * voxelSize);
  float3 camera_voxel = K * pos;
  if (pos.z < 0.0001f) // some near plane constraint
    return;
  // Interpolation version
  const float2 pixel = make_float2(camera_voxel.x / camera_voxel.z + 0.5f,
      camera_voxel.y / camera_voxel.z + 0.5f);
  if (pixel.x < 0.5f || pixel.x > depthSize.x - 1.5f || 
      pixel.y < 0.5f || pixel.y > depthSize.y - 1.5f) return;
  // const float depthSample = interpDepth(depth, depthSize, pixel);
  const uint2 px = make_uint2(pixel.x, pixel.y);
  const float depthSample = depth[px.x + depthSize.x*px.y];
  if (depthSample <=  0) return;

  const float diff = (camera_voxel.z - depthSample)
    * std::sqrt( 1 + sq(pos.x / pos.z) + sq(pos.y / pos.z));
  const float sample = HNew(diff/(noiseFactor*sq(camera_voxel.z)), camera_voxel.z);
  if(sample == 0.5f) return;
  typename VoxelBlock<BFusion>::compute_type data = node->value_; // should do an offsetted access here
  data.x = applyWindow(data.x, SURF_BOUNDARY, DELTA_T, CAPITAL_T);
  data.x = clamp(updateLogs(data.x, sample), BOTTOM_CLAMP, TOP_CLAMP);
  uint3 coords = unpack_morton(node->code);
  if((coords.x == 0) && coords.y == 256 && (coords.z == 0)) {
    printf("Vox: %d, %d, %d\n", x, y, z);
  }
  node ->value_ = data;
}

#endif
