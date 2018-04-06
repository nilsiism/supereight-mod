#ifndef MORTON_UTILS_HPP
#define MORTON_UTILS_HPP
#include <math_utils.h>
#include <cstdint>
#include <octree_defines.h>

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

inline uint3 unpack_morton(uint64_t code){
  return make_uint3(compact(code >> 0ull), compact(code >> 1ull), 
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

/*
 * Algorithm 5 of p4est paper: https://epubs.siam.org/doi/abs/10.1137/100791634
 */
inline uint3 face_neighbour(const octlib::key_t o, 
    const unsigned int face, const unsigned int l, 
    const unsigned int max_depth) {
  uint3 coords = unpack_morton(o);
  const unsigned int side = 1 << (max_depth - l); 
  coords.x = coords.x + ((face == 0) ? -side : (face == 1) ? side : 0);
  coords.y = coords.y + ((face == 2) ? -side : (face == 3) ? side : 0);
  coords.z = coords.z + ((face == 4) ? -side : (face == 5) ? side : 0);
  return {coords.x, coords.y, coords.z};
}

static inline void compute_prefix(const octlib::key_t * in, octlib::key_t * out,
    unsigned int num_keys, const octlib::key_t mask){

#pragma omp parallel for
  for (unsigned int i = 0; i < num_keys; i++){
    out[i] = in[i] & mask;
  }
}
#endif
