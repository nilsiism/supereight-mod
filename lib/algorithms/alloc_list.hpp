#ifndef ALLOC_LIST_HPP
#define ALLOC_LIST_HPP

#include <math_utils.h> 
#include <node.hpp>
#include <utils/morton_utils.hpp>

/* 
 * \brief Given a depth map and camera matrix it computes the list of 
 * voxels intersected but not allocated by the rays around the measurement m in
 * a region comprised between m +/- band. 
 * \param allocationList output list of keys corresponding to voxel blocks to
 * be allocated
 * \param reserved allocated size of allocationList
 * \param map_index indexing structure used to index voxel blocks 
 * \param pose camera extrinsics matrix
 * \param K camera intrinsics matrix
 * \param depthmap input depth map
 * \param imageSize dimensions of depthmap
 * \param size discrete extent of the map, in number of voxels
 * \param voxelSize spacing between two consegutive voxels, in metric space
 * \param band maximum extent of the allocating region, per ray
 */
template <typename FieldType, template <typename> class IndexType, typename HashType>
unsigned int buildAllocationList(HashType * allocationList, size_t reserved,
    IndexType<FieldType>& map_index, const Matrix4 &pose, const Matrix4& K, 
    const float *depthmap, const uint2 &imageSize, const unsigned int size,  
    const float voxelSize, const float band) {

  const float inverseVoxelSize = 1/voxelSize;

  Matrix4 invK = inverse(K);
  const Matrix4 kPose = pose * invK;

#ifdef _OPENMP
  std::atomic<unsigned int> voxelCount;
#else
  unsigned int voxelCount;
#endif

  unsigned int x, y;
  const float3 camera = get_translation(pose);
  const int numSteps = ceil(2*band*inverseVoxelSize);
  voxelCount = 0;
#pragma omp parallel for \
  private(y)
  for (y = 0; y < imageSize.y; y++) {
    for (x = 0; x < imageSize.x; x++) {
      if(depthmap[x + y*imageSize.x] == 0)
        continue;
      const float depth = depthmap[x + y*imageSize.x];
      float3 worldVertex = (kPose * make_float3((x + 0.5f) * depth, 
            (y + 0.5f) * depth, depth));

      float3 direction = normalize(worldVertex - camera);
      const float3 origin = worldVertex - band * direction;
      const float3 step = (2*direction*band)/numSteps;

      int3 voxel;
      float3 voxelPos = origin;
      for(int i = 0; i < numSteps; i++){
        float3 voxelScaled = floorf(voxelPos * inverseVoxelSize);
        if((voxelScaled.x < size) && (voxelScaled.y < size) &&
            (voxelScaled.z < size) && (voxelScaled.x >= 0) &&
            (voxelScaled.y >= 0) && (voxelScaled.z >= 0)){
          voxel = make_int3(voxelScaled);
          VoxelBlock<FieldType> * n = map_index.fetch(voxel.x, voxel.y, voxel.z);
          if(!n){
            HashType k = map_index.hash(voxel.x, voxel.y, voxel.z);
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
