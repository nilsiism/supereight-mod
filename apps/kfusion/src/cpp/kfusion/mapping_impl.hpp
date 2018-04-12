#ifndef KFUSION_MAPPING_HPP
#define KFUSION_MAPPING_HPP
#include <node.hpp>
#include "../continuous/volume_traits.hpp"

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
