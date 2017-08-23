#ifndef MORTON_UTILS_HPP
#define MORTON_UTILS_HPP
#include <math_utils.h>
#include <cstdint>

inline uint64_t expand(unsigned long long value) {
  uint64_t x = value & 0x1fffff;
  x = (x | x << 32) & 0x1f00000000ffff;
  x = (x | x << 16) & 0x1f0000ff0000ff;
  x = (x | x << 8)  & 0x100f00f00f00f00f;
  x = (x | x << 4)  & 0x10c30c30c30c30c3;
  x = (x | x << 2)  & 0x1249249249249249;
  return x;
}

inline uint64_t compact(uint64_t value) {
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

inline uint3 unpack_morton(uint64_t code){
  return make_uint3(compact(code >> 0), compact(code >> 1), 
                    compact(code >> 2));
}

inline unsigned long long compute_morton(uint64_t x, 
    uint64_t y, uint64_t z){
  uint64_t code = 0;

  x = expand(x);
  y = expand(y) << 1;
  z = expand(z) << 2;

  code = x | y | z;
  return code;
}

#endif
