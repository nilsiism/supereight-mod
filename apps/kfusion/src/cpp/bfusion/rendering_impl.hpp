#include <math_utils.h>
#include <type_traits>

float4 raycast(const Volume<BFusion>& volume, const uint2 pos, const Matrix4 view,
    const float nearPlane, const float farPlane, const float, 
    const float step, const float)
{
  const float3 origin = get_translation(view);
  const float3 direction = normalize(rotate(view, make_float3(pos.x, pos.y, 1.f)));

  // intersect ray with a box
  // http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm
  // compute intersection of ray with all six bbox planes
  const float3 invR = make_float3(1.0f) / direction;
  const float3 tbot = -1 * invR * origin;
  const float3 ttop = invR * (volume._dim - origin);

  // re-order intersections to find smallest and largest on each axis
  const float3 tmin = fminf(ttop, tbot);
  const float3 tmax = fmaxf(ttop, tbot);

  // find the largest tmin and the smallest tmax
  const float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y),
      fmaxf(tmin.x, tmin.z));
  const float smallest_tmax = fminf(fminf(tmax.x, tmax.y),
      fminf(tmax.x, tmax.z));

  // check against near and far plane
  const float tnear = fmaxf(largest_tmin, nearPlane);
  const float tfar = fminf(smallest_tmax, farPlane);
  auto select_occupancy = [](const auto& val){ return val.x; };

  if (tnear < tfar) {
    float t = tnear;
    float stepsize = step;
    float f_t = volume.interp(origin + direction * t, select_occupancy);
    float f_tt = 0;

    // if we are not already in it
    if (f_t <= SURF_BOUNDARY) { 
      for (; t < tfar; t += stepsize) {
        Volume<BFusion>::compute_type data = volume.get(origin + direction * t);
        if(data.y == 0){
          t += stepsize * 6.0f;
          continue;
        }
        if(data.x > -100.f && data.y > 0.f){
          f_tt = volume.interp(origin + direction * t, select_occupancy);
        }
        if (f_tt > SURF_BOUNDARY) break;
        f_t = f_tt;
      }            
      if (f_tt > SURF_BOUNDARY) {
        // got it, calculate accurate intersection
        t = t - stepsize * (f_tt - SURF_BOUNDARY) / (f_tt - f_t);
        return make_float4(origin + direction * t, t);
      }
    }
  }
  return make_float4(0);
}
