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

#ifndef MESHING_HPP
#define MESHING_HPP
#include <octree.hpp>
#include "edge_tables.h"

template <typename T>
float3 compute_intersection(Octree<T>& volume, const uint3 source, const uint3 dest){
 
  const float voxelSize = volume.dim()/volume.size(); 
  float3 s = make_float3(source.x * voxelSize, source.y * voxelSize, source.z * voxelSize);
  float3 d = make_float3(dest.x * voxelSize, dest.y * voxelSize, dest.z * voxelSize);
  float v1 = volume.get(source.x, source.y, source.z).x;
  float v2 = volume.get(dest.x, dest.y, dest.z).x; 
  return s + (0 - v1)*(d - s)/(v2-v1);
  
}

template <typename T>
float3 interp_vertexes(Octree<T>& volume, const uint x, const uint y, const uint z,
                          const int edge){
  switch(edge){
      case 0:  return compute_intersection(volume, make_uint3(x,   y, z),     
                                           make_uint3(x+1, y, z));
      case 1:  return compute_intersection(volume, make_uint3(x+1, y, z),     
                                           make_uint3(x+1, y, z+1));
      case 2:  return compute_intersection(volume, make_uint3(x+1, y, z+1),   
                                           make_uint3(x, y, z+1));
      case 3:  return compute_intersection(volume, make_uint3(x,   y, z),     
                                           make_uint3(x, y, z+1));
      case 4:  return compute_intersection(volume, make_uint3(x,   y+1, z),   
                                           make_uint3(x+1, y+1, z));
      case 5:  return compute_intersection(volume, make_uint3(x+1, y+1, z),   
                                           make_uint3(x+1, y+1, z+1));
      case 6:  return compute_intersection(volume, make_uint3(x+1, y+1, z+1), 
                                           make_uint3(x, y+1, z+1));
      case 7:  return compute_intersection(volume, make_uint3(x,   y+1, z),   
                                           make_uint3(x,   y+1, z+1));

      case 8:  return compute_intersection(volume, make_uint3(x,   y, z),     
                                           make_uint3(x,   y+1, z));
      case 9:  return compute_intersection(volume, make_uint3(x+1, y, z),     
                                           make_uint3(x+1, y+1, z));
      case 10: return compute_intersection(volume, make_uint3(x+1, y, z+1),   
                                           make_uint3(x+1, y+1, z+1));
      case 11: return compute_intersection(volume, make_uint3(x,   y, z+1),   
                                           make_uint3(x,   y+1, z+1));
      }
  return make_float3(0);
}

template <typename T>
uint8_t compute_index(const uint x, const uint y, const uint z,
                              VoxelBlock<T> *leaf){
   uint8_t index = 0;
   if(leaf->data(make_int3(x, y, z)).y == 0) 
     return 0;
   if(leaf->data(make_int3(x+1, y, z)).y == 0)
     return 0;
   if(leaf->data(make_int3(x+1, y, z+1)).y == 0)
     return 0;
   if(leaf->data(make_int3(x, y, z+1)).y == 0)
     return 0;
   if(leaf->data(make_int3(x, y+1, z)).y == 0)
     return 0;
   if(leaf->data(make_int3(x+1, y+1, z)).y == 0)
     return 0;
   if(leaf->data(make_int3(x+1, y+1, z+1)).y == 0)
     return 0;
   if(leaf->data(make_int3(x, y+1, z+1)).y == 0)
     return 0;

   if(leaf->data(make_int3(x, y, z)).x > 0.f) 
       index |= 1;
   if(leaf->data(make_int3(x+1, y, z)).x > 0.f)
       index |= 2;
   if(leaf->data(make_int3(x+1, y, z+1)).x > 0.f)
       index |= 4;
   if(leaf->data(make_int3(x, y, z+1)).x > 0.f)
       index |= 8;
   if(leaf->data(make_int3(x, y+1, z)).x > 0.f)
       index |= 16;
   if(leaf->data(make_int3(x+1, y+1, z)).x > 0.f)
       index |= 32;
   if(leaf->data(make_int3(x+1, y+1, z+1)).x > 0.f)
       index |= 64;
   if(leaf->data(make_int3(x, y+1, z+1)).x > 0.f)
       index |= 128;
   return index;
}

template <typename T>
uint8_t compute_index_at_border(Octree<T>& volume, const uint x, const uint y, 
                                           const uint z){
   uint8_t index = 0;
   // constexpr auto& volume.get = Octree<T>::volume.get;
   if(volume.get(x, y, z).y == 0) 
     return 0;
   if(volume.get(x+1, y, z).y == 0)
     return 0;
   if(volume.get(x+1, y, z+1).y == 0)
     return 0;
   if(volume.get(x, y, z+1).y == 0)
     return 0;
   if(volume.get(x, y+1, z).y == 0)
     return 0;
   if(volume.get(x+1, y+1, z).y == 0)
     return 0;
   if(volume.get(x+1, y+1, z+1).y == 0)
     return 0;
   if(volume.get(x, y+1, z+1).y == 0)
     return 0;

   if(volume.get(x, y, z).x > 0.f) 
       index |= 1;
   if(volume.get(x+1, y, z).x > 0.f)
       index |= 2;
   if(volume.get(x+1, y, z+1).x > 0.f)
       index |= 4;
   if(volume.get(x, y, z+1).x > 0.f)
       index |= 8;
   if(volume.get(x, y+1, z).x > 0.f)
       index |= 16;
   if(volume.get(x+1, y+1, z).x > 0.f)
       index |= 32;
   if(volume.get(x+1, y+1, z+1).x > 0.f)
       index |= 64;
   if(volume.get(x, y+1, z+1).x > 0.f)
       index |= 128;
   return index;
}

namespace algorithms {
  template <typename T, typename TriangleType>
    void marching_cube(Octree<T>& volume, std::vector<TriangleType>& triangles)
    {

      std::stringstream points, polygons;
      std::vector<VoxelBlock<T>*> blocklist;
      std::mutex lck;
      const int size = volume.size();
      volume.getBlockList(blocklist, false);
      std::cout << "Blocklist size: " << blocklist.size() << std::endl;

#pragma omp parallel for
      for(size_t i = 0; i < blocklist.size(); i++){
        VoxelBlock<T> * leaf = static_cast<VoxelBlock<T> *>(blocklist[i]);  
        int edge = VoxelBlock<T>::side;
        int x, y, z ; 
        int xbound = clamp(leaf->coordinates().x + edge, 0, size-1);
        int ybound = clamp(leaf->coordinates().y + edge, 0, size-1);
        int zbound = clamp(leaf->coordinates().z + edge, 0, size-1);
        for(x = leaf->coordinates().x; x < xbound; x++){
          for(y = leaf->coordinates().y; y < ybound; y++){
            for(z = leaf->coordinates().z; z < zbound; z++){

              uint8_t index = (x == xbound-1 || y == ybound-1 || z == zbound-1 ) ?
                compute_index_at_border(volume, x, y, z) :
                compute_index(x, y, z, leaf);

              int * edges = triTable[index]; 
              for(unsigned int e = 0; edges[e] != -1 && e < 16; e += 3){
                float3 v1 = interp_vertexes(volume, x, y, z, edges[e]);
                float3 v2 = interp_vertexes(volume, x, y, z, edges[e+1]);
                float3 v3 = interp_vertexes(volume, x, y, z, edges[e+2]);
                Triangle temp = Triangle();
                temp.vertexes[0] = v1;
                temp.vertexes[1] = v2;
                temp.vertexes[2] = v3;
                lck.lock();
                triangles.push_back(temp);
                lck.unlock();
              }
            }
          }
        }
      }
    }
}
#endif
