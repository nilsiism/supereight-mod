#ifndef KFUSION_MAPPING_HPP
#define KFUSION_MAPPING_HPP
#include <se/node.hpp>

struct sdf_update {

  template <typename DataHandlerT>
  void operator()(DataHandlerT& handler, const Eigen::Vector3i&, 
      const Eigen::Vector3f& pos, const Eigen::Vector2f& pixel) {

    const Eigen::Vector2i px = pixel.cast<int>();
    const float depthSample = depth[px(0) + depthSize(0)*px(1)];
    if (depthSample <=  0) return;
    const float diff = (depthSample - pos(2)) 
      * std::sqrt( 1 + sq(pos(0) / pos(2)) + sq(pos(1) / pos(2)));
    if (diff > -mu) {
      const float sdf = fminf(1.f, diff / mu);
      auto data = handler.get();
      data.x = clamp((data.y * data.x + sdf) / (data.y + 1), -1.f,
          1.f);
      data.y = fminf(data.y + 1, maxweight);
      handler.set(data);
    }
  } 

  sdf_update(const float * d, const Eigen::Vector2i framesize, float m, int mw) : 
    depth(d), depthSize(framesize), mu(m), maxweight(mw){};

  const float * depth;
  Eigen::Vector2i depthSize;
  float mu;
  int maxweight;
};

#endif
