#ifndef KFUSION_MAPPING_HPP
#define KFUSION_MAPPING_HPP
#include <node.hpp>
#include "../continuous/volume_traits.hpp"

void integrate(VoxelBlock<SDF> * block, const float * depth, uint2 depthSize, 
    const float voxelSize, const Matrix4 invTrack, const Matrix4 K, 
    const float mu, const float maxweight) { 

  using VoxelBlockT = VoxelBlock<SDF>;
  const float3 delta = rotate(invTrack, make_float3(voxelSize, 0, 0));
  const float3 cameraDelta = rotate(K, delta);
  const int3 blockCoord = block->coordinates();
  bool is_visible = false;
  block->active(is_visible);

  unsigned int y, z, blockSide; 
  blockSide = VoxelBlockT::side;
  unsigned int ylast = blockCoord.y + blockSide;
  unsigned int zlast = blockCoord.z + blockSide;
  for(z = blockCoord.z; z < zlast; ++z)
    for (y = blockCoord.y; y < ylast; ++y){
      int3 pix = make_int3(blockCoord.x, y, z);
      float3 start = invTrack * make_float3((pix.x + 0.5f) * voxelSize, 
          (pix.y + 0.5f) * voxelSize, (pix.z + 0.5f) * voxelSize);
			float3 camerastart = K * start;
#pragma simd
      for (unsigned int x = 0; x < blockSide; ++x){
        pix.x = x + blockCoord.x; 
        const float3 cameraX = camerastart + (x*cameraDelta);
        const float3 pos = start + (x*delta);
        if (pos.z < 0.0001f) // some near plane constraint
          continue;
        const float2 pixel = make_float2(cameraX.x / cameraX.z + 0.5f,
            cameraX.y / cameraX.z + 0.5f);
        if (pixel.x < 0 || pixel.x > depthSize.x - 1 || pixel.y < 0
            || pixel.y > depthSize.y - 1)
          continue;
        const uint2 px = make_uint2(pixel.x, pixel.y);
        if (depth[px.x + px.y * depthSize.x] == 0)
          continue;
        is_visible = true;
        const float diff =
          (depth[px.x + px.y * depthSize.x] - cameraX.z)
          * std::sqrt( 1 + sq(pos.x / pos.z) + sq(pos.y / pos.z));
        if (diff > -mu) {
          const float sdf = fminf(1.f, diff / mu);
          typename VoxelBlockT::compute_type data = block->data(pix); // should do an offsetted access here
          data.x = clamp((data.y * data.x + sdf) / (data.y + 1), -1.f,
              1.f);
          data.y = fminf(data.y + 1, maxweight);
          block->data(pix, data);
        }
      }
    }
  block->active(is_visible);
}

struct sdf_update {

  template <typename DataHandlerT>
  void operator()(DataHandlerT& handler, const int3&, const float3& pos, 
     const float2& pixel) {

    const uint2 px = make_uint2(pixel.x, pixel.y);
    const float depthSample = depth[px.x + depthSize.x*px.y];
    if (depthSample <=  0) return;
    const float diff = (depthSample - pos.z) 
      * std::sqrt( 1 + sq(pos.x / pos.z) + sq(pos.y / pos.z));
    if (diff > -mu) {
      const float sdf = fminf(1.f, diff / mu);
      auto data = handler.get();
      data.x = clamp((data.y * data.x + sdf) / (data.y + 1), -1.f,
          1.f);
      data.y = fminf(data.y + 1, maxweight);
      handler.set(data);
    }
  } 

  sdf_update(const float * d, const uint2 framesize, float m, int mw) : 
    depth(d), depthSize(framesize), mu(m), maxweight(mw){};

  const float * depth;
  uint2 depthSize;
  float mu;
  int maxweight;
};

#endif
