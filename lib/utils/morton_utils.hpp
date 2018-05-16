#ifndef MORTON_UTILS_HPP
#define MORTON_UTILS_HPP
#include <cstdint>
#include "octree_defines.h"
#include "se_common.h"

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

inline Eigen::Vector3i unpack_morton(uint64_t code){
  return Eigen::Vector3i(compact(code >> 0ull), compact(code >> 1ull), 
                    compact(code >> 2ull));
}

inline uint64_t compute_morton(uint64_t x, 
    uint64_t y, uint64_t z){
  uint64_t code = 0;

  x = expand(x);
  y = expand(y) << 1;
  z = expand(z) << 2;

  code = x | y | z;
  return code;
}

static inline void compute_prefix(const octlib::key_t * in, octlib::key_t * out,
    unsigned int num_keys, const octlib::key_t mask){

#pragma omp parallel for
  for (unsigned int i = 0; i < num_keys; i++){
    out[i] = in[i] & mask;
  }
}
#endif
