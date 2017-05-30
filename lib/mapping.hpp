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

#ifndef MAPPING_HPP
#define MAPPING_HPP
#include <node.hpp>

float3 voxelToPos(const uint3 p, const float voxelSize){
  return make_float3((p.x + 0.5f) * voxelSize, (p.y + 0.5f) * voxelSize,
    (p.z + 0.5f) * voxelSize);
}
template <typename T>
void integratePass(Aggregate<T> ** blockList, unsigned int list_size, 
    const float * depth, uint2 depthSize, const float voxelSize, 
    const Matrix4 invTrack, const Matrix4 K, const float mu, 
    const float maxweight, const int current_frame) {

#pragma omp parallel for
  for(int i = 0; i < list_size; ++i){
      integrate(blockList[i], depth, depthSize, voxelSize,
          invTrack, K, mu, maxweight);
      blockList[i]->last_integrated_frame(current_frame);
  }
}

template <typename T>
void integrate(Aggregate<T> * block, const float * depth, uint2 depthSize, 
    const float voxelSize, const Matrix4 invTrack, const Matrix4 K, 
    const float mu, const float maxweight) { 

  const float3 delta = rotate(invTrack, make_float3(0, 0, voxelSize));
  const float3 cameraDelta = rotate(K, delta);
  const uint3 blockCoord = block->coordinates();
  bool is_visible = false;
  block->active(is_visible);

  unsigned int x, y, z, blockSide; 
  blockSide = Aggregate<T>::side;
  unsigned int xlast = blockCoord.x + blockSide;
  unsigned int ylast = blockCoord.y + blockSide;
  unsigned int zlast = blockCoord.z + blockSide;
  for (y = blockCoord.y; y < ylast; ++y){
    for (x = blockCoord.x; x < xlast; ++x){
      uint3 pix = make_uint3(x, y, blockCoord.z);
      float3 pos = invTrack * voxelToPos(pix, voxelSize);
      float3 cameraX = K * pos;

      for (; pix.z < zlast; ++pix.z, pos+= delta, cameraX += cameraDelta){
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
          typename Aggregate<T>::compute_type data = block->data(pix); // should do an offsetted access here
          data.x = clamp((data.y * data.x + sdf) / (data.y + 1), -1.f,
              1.f);
          data.y = fminf(data.y + 1, maxweight);
          block->data(pix, data);
        }
      }
    }
  }
  block->active(is_visible);
}

#endif
