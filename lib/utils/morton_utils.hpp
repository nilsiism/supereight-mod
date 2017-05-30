#ifndef MORTON_UTILS_HPP
#define MORTON_UTILS_HPP
#include <math_utils.h>

inline unsigned int expand(uint value){
  unsigned int x;
  x = value & 0x03FF;
  x = ((x << 16) + x) & 0xFF0000FF;
  x = ((x <<  8) + x) & 0x0F00F00F;
  x = ((x <<  4) + x) & 0xC30C30C3;
  x = ((x <<  2) + x) & 0x49249249;
  return x;
}

inline uint compact(uint value){
    unsigned int x;
    x = value & 0x49249249;
    x = (x ^ (x >>  2)) & 0xC30C30C3;
    x = (x ^ (x >>  4)) & 0x0F00F00F;
    x = (x ^ (x >>  8)) & 0xFF0000FF;
    x = (x ^ (x >> 16)) & 0x03FF;
    return x;
}

inline uint3 unpack_morton(uint code){
  return make_uint3(compact(code >> 0), compact(code >> 1), 
                    compact(code >> 2));
}

inline uint compute_morton(const uint2 coordinates){
  uint x, y = 0;
  uint code = 0;

  x = expand(coordinates.x);
  y = expand(coordinates.y);

  code = x | (y << 1);
  return code;
}

inline uint compute_morton(const uint3 coordinates){
  uint x, y, z = 0;
  uint code = 0;

  x = expand(coordinates.x);
  y = expand(coordinates.y);
  z = expand(coordinates.z);

  code = x | (y << 1 | (z << 2));
  return code;
}

inline uint compute_morton(const int3 coordinates){
  uint x, y, z = 0;
  uint code = 0;

  x = expand(coordinates.x);
  y = expand(coordinates.y);
  z = expand(coordinates.z);

  code = x | (y << 1 | (z << 2));
  return code;
}

inline uint compute_morton(const uint xi, const uint yi, const uint zi){
  uint x, y, z = 0;
  uint code = 0;

  x = expand(xi);
  y = expand(yi);
  z = expand(zi);

  code = x | (y << 1 | (z << 2));
  return code;
}

inline int compute_morton(const int xi, const int yi, const int zi){
  int x, y, z = 0;
  int code = 0;

  x = expand(xi);
  y = expand(yi);
  z = expand(zi);

  code = x | (y << 1 | (z << 2));
  return code;
}
#endif
