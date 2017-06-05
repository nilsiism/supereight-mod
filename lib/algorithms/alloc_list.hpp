#ifndef ALLOC_LIST_HPP
#define ALLOC_LIST_HPP

#include <math_utils.h> 
#include <node.hpp>
#include <algorithms/mapping.hpp>
#include <utils/morton_utils.hpp>

template <typename FieldType, template <typename> class IndexType>
unsigned int buildAllocationList(uint * allocationList, size_t reserved,
    IndexType<FieldType>& map_index, const Matrix4 &pose, const Matrix4& K, 
    const float *depthmap, const uint2 &imageSize, const unsigned int size,  
    const float voxelSize, const float mu, const int frame) {

  const float inverseVoxelSize = 1/voxelSize;

  const Matrix4 invPose = inverse(pose); Matrix4 invK = inverse(K);
  const Matrix4 kPose = pose * invK;

#ifdef _OPENMP
  std::atomic<unsigned int> voxelCount;
#else
  unsigned int voxelCount;
#endif

  unsigned int x, y;
  const float3 camera = get_translation(pose);
  const int numSteps = ceil(2*mu*inverseVoxelSize);
  voxelCount = 0;
#pragma omp parallel for \
  private(y)
  for (y = 0; y < imageSize.y; y++) {
    for (x = 0; x < imageSize.x; x++) {
      if(depthmap[x + y*imageSize.x] == 0)
        continue;
      const float depth = depthmap[x + y*imageSize.x];
      float3 worldVertex = (kPose * make_float3(x * depth, y * depth, depth));

      float3 direction = normalize(worldVertex - camera);
      const float3 origin = worldVertex - mu * direction;
      const float3 end = worldVertex + mu * direction;
      const float3 step = (2*direction*mu)/numSteps;

      int3 voxel;
      float3 voxelPos = origin;
      for(int i = 0; i < numSteps; i++){
        float3 voxelScaled = floorf(voxelPos * inverseVoxelSize);
        if((voxelScaled.x < size) && (voxelScaled.y < size) &&
            (voxelScaled.z < size) && (voxelScaled.x >= 0) &&
            (voxelScaled.y >= 0) && (voxelScaled.z >= 0)){
          voxel = make_int3(voxelScaled);
          Aggregate<FieldType> * n = map_index.fetch(voxel.x, voxel.y, voxel.z);
          if(!n){
            uint k = map_index.hash(voxel.x, voxel.y, voxel.z);
            unsigned int idx = ++voxelCount;
            if(idx < reserved) {
              allocationList[idx] = k;
            } else
              break;
          }
          else {
            n->active(true); 
          }
        }
        voxelPos +=step;
      }
    }
  }

  const uint written = voxelCount;
  return written >= reserved ? reserved : written;
}
#endif
