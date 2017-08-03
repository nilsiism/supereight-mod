#ifndef MORTON_UTILS_HPP
#define MORTON_UTILS_HPP
#include <math_utils.h>

// inline unsigned int expand(unsigned int value){
//   unsigned int x;
//   x = value & 0x03FF;
//   x = ((x << 16) + x) & 0xFF0000FF;
//   x = ((x <<  8) + x) & 0x0F00F00F;
//   x = ((x <<  4) + x) & 0xC30C30C3;
//   x = ((x <<  2) + x) & 0x49249249;
//   return x;
// }

inline unsigned long long expand(unsigned long long value) {
  uint64_t x = value & 0x1fffff;
  x = (x | x << 32) & 0x1f00000000ffff;
  x = (x | x << 16) & 0x1f0000ff0000ff;
  x = (x | x << 8)  & 0x100f00f00f00f00f;
  x = (x | x << 4)  & 0x10c30c30c30c30c3;
  x = (x | x << 2)  & 0x1249249249249249;
  return x;
}

inline unsigned long long compact(unsigned long long value) {
  uint64_t x = value & 0x1249249249249249;
  x = (x | x >> 2)   & 0x10c30c30c30c30c3;
  x = (x | x >> 4)   & 0x100f00f00f00f00f;
  x = (x | x >> 8)   & 0x1f0000ff0000ff;
  x = (x | x >> 16)  & 0x1f00000000ffff;
  x = (x | x >> 32)  & 0x1fffff;
  return x;
}

// inline uint compact(unsigned int value){
//     unsigned int x;
//     x = value & 0x49249249;
//     x = (x ^ (x >>  2)) & 0xC30C30C3;
//     x = (x ^ (x >>  4)) & 0x0F00F00F;
//     x = (x ^ (x >>  8)) & 0xFF0000FF;
//     x = (x ^ (x >> 16)) & 0x03FF;
//     return x;
// }

inline uint3 unpack_morton(unsigned long long code){
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

inline unsigned long long compute_morton(unsigned int x, int y,
    int z){
  unsigned long long int code = 0;

  x = expand(x);
  y = expand(y);
  z = expand(z);

  code = x | (y << 1 | (z << 2));
  return code;
}

#endif
