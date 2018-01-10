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

namespace meshing {
  enum status : uint8_t {
    OUTSIDE = 0x0,
    UNKNOWN = 0xFE, // 254
    INSIDE = 0xFF, // 255
  };

  template <typename Map, typename FieldSelector>
    inline float3 compute_intersection(const Map& volume, FieldSelector select,
        const uint3& source, const uint3& dest){
      const float voxelSize = volume.dim()/volume.size(); 
      float3 s = make_float3(source.x * voxelSize, source.y * voxelSize, source.z * voxelSize);
      float3 d = make_float3(dest.x * voxelSize, dest.y * voxelSize, dest.z * voxelSize);
      float v1 = select(volume.get_fine(source.x, source.y, source.z));
      float v2 = select(volume.get_fine(dest.x, dest.y, dest.z)); 
      return s + (0.0 - v1)*(d - s)/(v2-v1);
    }

  template <typename Map, typename FieldSelector>
    inline float3 interp_vertexes(const Map& volume, FieldSelector select, 
        const uint x, const uint y, const uint z, const int edge){
      switch(edge){
        case 0:  return compute_intersection(volume, select, make_uint3(x,   y, z),     
                     make_uint3(x+1, y, z));
        case 1:  return compute_intersection(volume, select, make_uint3(x+1, y, z),     
                     make_uint3(x+1, y, z+1));
        case 2:  return compute_intersection(volume, select, make_uint3(x+1, y, z+1),   
                     make_uint3(x, y, z+1));
        case 3:  return compute_intersection(volume, select, make_uint3(x,   y, z),     
                     make_uint3(x, y, z+1));
        case 4:  return compute_intersection(volume, select, make_uint3(x,   y+1, z),   
                     make_uint3(x+1, y+1, z));
        case 5:  return compute_intersection(volume, select, make_uint3(x+1, y+1, z),   
                     make_uint3(x+1, y+1, z+1));
        case 6:  return compute_intersection(volume, select, make_uint3(x+1, y+1, z+1), 
                     make_uint3(x, y+1, z+1));
        case 7:  return compute_intersection(volume, select, make_uint3(x,   y+1, z),   
                     make_uint3(x,   y+1, z+1));

        case 8:  return compute_intersection(volume, select, make_uint3(x,   y, z),     
                     make_uint3(x,   y+1, z));
        case 9:  return compute_intersection(volume, select, make_uint3(x+1, y, z),     
                     make_uint3(x+1, y+1, z));
        case 10: return compute_intersection(volume, select, make_uint3(x+1, y, z+1),   
                     make_uint3(x+1, y+1, z+1));
        case 11: return compute_intersection(volume, select, make_uint3(x,   y, z+1),   
                     make_uint3(x,   y+1, z+1));
      }
      return make_float3(0);
    }

  template <typename FieldType, typename PointT>
    inline void gather_points( const VoxelBlock<FieldType>* cached, PointT points[8], 
        const int x, const int y, const int z) {
      points[0] = cached->data(make_int3(x, y, z)); 
      points[1] = cached->data(make_int3(x+1, y, z));
      points[2] = cached->data(make_int3(x+1, y, z+1));
      points[3] = cached->data(make_int3(x, y, z+1));
      points[4] = cached->data(make_int3(x, y+1, z));
      points[5] = cached->data(make_int3(x+1, y+1, z));
      points[6] = cached->data(make_int3(x+1, y+1, z+1));
      points[7] = cached->data(make_int3(x, y+1, z+1));
    }

  template <typename FieldType, template <typename FieldT> class MapT, typename PointT>
  inline void gather_points(const MapT<FieldType>& volume, PointT points[8], 
                 const int x, const int y, const int z) {
               points[0] = volume.get_fine(x, y, z); 
               points[1] = volume.get_fine(x+1, y, z);
               points[2] = volume.get_fine(x+1, y, z+1);
               points[3] = volume.get_fine(x, y, z+1);
               points[4] = volume.get_fine(x, y+1, z);
               points[5] = volume.get_fine(x+1, y+1, z);
               points[6] = volume.get_fine(x+1, y+1, z+1);
               points[7] = volume.get_fine(x, y+1, z+1);
             }

  template <typename FieldType, template <typename FieldT> class MapT,
  typename InsidePredicate>
  uint8_t compute_index(const MapT<FieldType>& volume, 
  const VoxelBlock<FieldType>* cached, InsidePredicate inside,
  const uint x, const uint y, const uint z){
    unsigned int blockSize =  VoxelBlock<FieldType>::side;
    unsigned int local = ((x % blockSize == blockSize - 1) << 2) | 
      ((y % blockSize == blockSize - 1) << 1) |
      ((z % blockSize) == blockSize - 1);

    typename MapT<FieldType>::compute_type points[8];
    if(!local) gather_points(cached, points, x, y, z);
    else gather_points(volume, points, x, y, z);

    uint8_t index = 0;

    // if(points[0].y == 0.f) return 0;
    // if(points[1].y == 0.f) return 0;
    // if(points[2].y == 0.f) return 0;
    // if(points[3].y == 0.f) return 0;
    // if(points[4].y == 0.f) return 0;
    // if(points[5].y == 0.f) return 0;
    // if(points[6].y == 0.f) return 0;
    // if(points[7].y == 0.f) return 0;

    std::cerr << points[0].x << " " << points[1].x << " " << points[2].x << " " 
              << points[3].x << " " << points[4].x << " " << points[5].x << " " 
              <<  points[6].x << " " << points[7].x << std::endl;

    if(inside(points[0])) index |= 1;
    if(inside(points[1])) index |= 2;
    if(inside(points[2])) index |= 4;
    if(inside(points[3])) index |= 8;
    if(inside(points[4])) index |= 16;
    if(inside(points[5])) index |= 32;
    if(inside(points[6])) index |= 64;
    if(inside(points[7])) index |= 128;
    std::cerr << std::endl << std::endl;

    return index;
  }

  inline bool checkVertex(const float3 v, const int dim){
    return (v.x <= 0 || v.y <=0 || v.z <= 0 || v.x > dim || v.y > dim || v.z > dim);
  }

}
namespace algorithms {
  template <typename FieldType, typename FieldSelector, typename TriangleType>
    void marching_cube(Octree<FieldType>& volume, FieldSelector select,
        std::vector<TriangleType>& triangles)
    {

      using namespace meshing;
      std::stringstream points, polygons;
      std::vector<VoxelBlock<FieldType>*> blocklist;
      std::mutex lck;
      const int size = volume.size();
      const int dim = volume.dim();
      volume.getBlockList(blocklist, false);
      std::cout << "Blocklist size: " << blocklist.size() << std::endl;
      

#pragma omp parallel for
      for(size_t i = 0; i < blocklist.size(); i++){
        VoxelBlock<FieldType> * leaf = static_cast<VoxelBlock<FieldType> *>(blocklist[i]);  
        int edge = VoxelBlock<FieldType>::side;
        int x, y, z ; 
        int xbound = clamp(leaf->coordinates().x + edge, 0, size-1);
        int ybound = clamp(leaf->coordinates().y + edge, 0, size-1);
        int zbound = clamp(leaf->coordinates().z + edge, 0, size-1);
        for(x = leaf->coordinates().x; x < xbound; x++){
          for(y = leaf->coordinates().y; y < ybound; y++){
            for(z = leaf->coordinates().z; z < zbound; z++){

              uint8_t index = meshing::compute_index(volume, leaf, select, x, y, z);

                int * edges = triTable[index]; 
              for(unsigned int e = 0; edges[e] != -1 && e < 16; e += 3){
                float3 v1 = interp_vertexes(volume, select, x, y, z, edges[e]);
                float3 v2 = interp_vertexes(volume, select, x, y, z, edges[e+1]);
                float3 v3 = interp_vertexes(volume, select, x, y, z, edges[e+2]);
                if(checkVertex(v1, dim) || checkVertex(v2, dim) || checkVertex(v3, dim)) continue;
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
