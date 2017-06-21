#ifndef BFUSION_MAPPING_HPP
#define BFUSION_MAPPING_HPP

#include <node.hpp>
#include <constant_parameters.h>

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

  static constexpr float interp_thresh = 0.05f;
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


inline float3 voxelToPos(const uint3 p, const float voxelSize){
  return make_float3((p.x + 0.5f) * voxelSize, (p.y + 0.5f) * voxelSize,
    (p.z + 0.5f) * voxelSize);
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
constexpr float scale_factor = (1.f - (farPlane - nearPlane) * const_offset); 

static inline float const_offset_integral(float t){
  float value = 0.f;
  if (nearPlane <= t && t <= farPlane)
      return (t - nearPlane) * const_offset;
  else if (farPlane < t)
      return (farPlane - nearPlane) * const_offset;
  return value;
}

static inline float HNew(const float val,const  float d_xr){
  const float Q_1 = bspline(val)    * scale_factor + const_offset_integral(d_xr      );
  const float Q_2 = bspline(val - 3)* scale_factor + const_offset_integral(d_xr - 3.f);
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
    const float noiseFactor, const float ) {

  const float3 delta = rotate(world2cam, make_float3(voxelSize, 0, 0));
  const float3 cameraDelta = rotate(K, delta);
  const uint3 blockCoord = block->coordinates();
  bool is_visible = false;
  block->active(is_visible);

  unsigned int y, z, blockSide; 
  blockSide = VoxelBlock<BFusion>::side;
  unsigned int ylast = blockCoord.y + blockSide;
  unsigned int zlast = blockCoord.z + blockSide;
  for(z = blockCoord.z; z < zlast; ++z)
    for (y = blockCoord.y; y < ylast; ++y){
      uint3 pix = make_uint3(blockCoord.x, y, z);
      float3 start = world2cam * make_float3((pix.x + 0.5f) * voxelSize, 
          (pix.y + 0.5f) * voxelSize, (pix.z + 0.5f) * voxelSize);
      float3 camerastart = K * start;
#pragma simd
      for (unsigned int x = 0; x < blockSide; ++x){
        pix.x = x + blockCoord.x; 
        const float3 camera_voxel = camerastart + (x*cameraDelta);
        const float3 pos = start + (x*delta);
        if (pos.z < 0.0001f) continue;

        // Original version
        // const float2 pixel = make_float2( camera_voxel.x / camera_voxel.z + 0.5f,
        //                                   camera_voxel.y / camera_voxel.z + 0.5f);
        // if (pixel.x < 0 || pixel.x > depth.size.x - 1 || 
        //     pixel.y < 0 || pixel.y > depth.size.y - 1) continue;
        // const uint2 px = make_uint2(pixel.x, pixel.y);
        // const float depthSample = depth[px];
        // if (depthSample ==  0) continue;

        // Interpolation version
        const float2 pixel = make_float2( camera_voxel.x / camera_voxel.z,
            camera_voxel.y / camera_voxel.z);
        if (pixel.x < 0.5f || pixel.x > depthSize.x - 1.5f || 
            pixel.y < 0.5f || pixel.y > depthSize.y - 1.5f) continue;
        const float depthSample = interpDepth(depth, depthSize, pixel);
        if (depthSample <=  0) continue;

        // const float diff = (camera_voxel.z - depthSample)/(noiseFactor*sq(camera_voxel.z));
        const float diff = (camera_voxel.z - depthSample)/(noiseFactor);
        const float sample = HNew(diff, camera_voxel.z);

        typename VoxelBlock<BFusion>::compute_type data = block->data(pix); // should do an offsetted access here
        data.x = clamp(updateLogs(data.x, sample), BOTTOM_CLAMP, TOP_CLAMP);
        data.x = applyWindow(data.x, SURF_BOUNDARY, DELTA_T, CAPITAL_T);
        block->data(pix, data);
      }
    }
}

#endif
