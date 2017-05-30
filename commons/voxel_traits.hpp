/*

Copyright 2016 Emanuele Vespa, Imperial College London 

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

*/

#ifndef _VOXEL_TRAITS_
#define _VOXEL_TRAITS_
#include <vector_types.h>
#include <cutil_math.h>
class SDF;
class SigmaSDF;
class ColouredSDF;

template <class VoxelTraits>
struct kfusion_voxel_traits{ };

template<>
struct kfusion_voxel_traits<SDF> {
  typedef float2 ComputeType;
  typedef short2 StoredType;
  static inline ComputeType empty(){ return make_float2(1.f, -1.f); }
  static inline StoredType initValue(){ return make_short2(32766, 0); }
  static inline StoredType update(const ComputeType value) {
     return make_short2(value.x * 32766.0f, value.y);
  }
  static inline ComputeType translate(const StoredType value){
    return make_float2(value.x * 0.00003051944088f, value.y);
  }
};

template<>
struct kfusion_voxel_traits<SigmaSDF> {
  typedef float ComputeType;
  typedef short StoredType;
  static inline ComputeType empty(){ return 0.5f; }
  static inline StoredType initValue(){ return 0.5f * 32766.0f ; }
  static inline StoredType update(const ComputeType value) {
     return value * 32766.0f;
  }
  static inline ComputeType translate(const StoredType value){
    return value * 0.00003051944088f;
  }
};

struct SDFVoxelRGB {
  float x;
  float y;
  unsigned char r;
  unsigned char g;
  unsigned char b;
  unsigned char w;
};

struct SDFVoxelRGBCompressed {
  short x;
  short y;
  uint8_t r;
  uint8_t g;
  uint8_t b;
  short w;
};

template<>
struct kfusion_voxel_traits<ColouredSDF> {
  typedef SDFVoxelRGB ComputeType ;
  typedef SDFVoxelRGBCompressed StoredType;
  static inline ComputeType empty(){ return {1.f, -1.f, 0, 0, 0}; }
  static inline StoredType initValue(){ return {1,0,0,0}; }
  static inline StoredType update(const ComputeType value) {
    StoredType temp;
    temp.x = value.x * 32766.0f;
    temp.y = value.y;
    temp.r = value.r;
    temp.g = value.g;
    temp.b = value.b;
    temp.w = value.w * 32766.f;
    return temp;
  }
  static inline ComputeType translate(const StoredType value){
    ComputeType temp;
    temp.x = value.x * 0.00003051944088f;
    temp.y = value.y;
    temp.r = value.r;
    temp.g = value.g;
    temp.b = value.b;
    temp.w = value.w * 0.00003051944088f;
    return temp;
  }
};
#endif
