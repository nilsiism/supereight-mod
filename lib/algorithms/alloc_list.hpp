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
template <typename FieldType, template <typename> class IndexType>
unsigned int buildAllocationList(uint * allocationList, size_t reserved,
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

template <typename FieldType, template <typename> class IndexType>
unsigned int buildIntersectionList(uint * allocationList, size_t reserved,
    IndexType<FieldType>& map_index, const Matrix4 &pose, const Matrix4& K, 
    const float *depthmap, const uint2 &imageSize, const unsigned int tree_depth,  
    const float voxelSize, const float stepsize, const float band) {

  const float inverseVoxelSize = 1.f/voxelSize;
  Matrix4 invK = inverse(K);
  const Matrix4 kPose = pose * invK;
  const int size = map_index.size();

#ifdef _OPENMP
  std::atomic<unsigned int> voxelCount;
#else
  unsigned int voxelCount;
#endif

  unsigned int x, y;
  const float3 camera = get_translation(pose);
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
      const float3 origin = camera;
      const float3 end = worldVertex - band * direction;
      const int numSteps = ceil(length(end - camera)/stepsize);
      const float3 step = direction*stepsize;

      int3 voxel;
      float3 voxelPos = origin;
      for(int i = 0; i < numSteps; i++){
        float3 voxelScaled = floorf(voxelPos * inverseVoxelSize);
        if((voxelScaled.x < size) && (voxelScaled.y < size) &&
            (voxelScaled.z < size) && (voxelScaled.x >= 0) &&
            (voxelScaled.y >= 0) && (voxelScaled.z >= 0)){
          voxel = make_int3(voxelScaled);
          if(!map_index.fetch_octant(voxel.x, voxel.y, voxel.z, tree_depth)){
            uint k = map_index.hash(voxel.x, voxel.y, voxel.z);
            unsigned int idx = ++voxelCount;
            if(idx < reserved) {
              allocationList[idx] = k;
            }
          }
        }
        voxelPos +=step;
      }
    }
  }
  const uint written = voxelCount;
  return written >= reserved ? reserved : written;
}

template <typename FieldType, template <typename> class IndexType>
void buildOctantList(uint * allocationList[2], size_t reserved[2],
    size_t written[2], IndexType<FieldType>& map_index, const Matrix4 &pose, 
    const Matrix4& K, const float *depthmap, const uint2 &imageSize, 
    const unsigned int tree_depth[2],  const float voxelSize, 
    const float stepsize[2], const float band) {

  const float inverseVoxelSize = 1.f/voxelSize;
  Matrix4 invK = inverse(K);
  const Matrix4 kPose = pose * invK;
  const int size = map_index.size();

#ifdef _OPENMP
  std::atomic<unsigned int> voxelCount[2];
#else
  unsigned int voxelCount[2];
#endif

  unsigned int x, y;
  const float3 camera = get_translation(pose);
  voxelCount[0] = 0;
  voxelCount[1] = 0;
#pragma omp parallel for \
  private(y)
  for (y = 0; y < imageSize.y; y++) {
    for (x = 0; x < imageSize.x; x++) {
      if(depthmap[x + y*imageSize.x] == 0)
        continue;
      int step_phase = 0;
      const float depth = depthmap[x + y*imageSize.x];
      float3 worldVertex = (kPose * make_float3((x + 0.5f) * depth, 
            (y + 0.5f) * depth, depth));

      float3 direction = normalize(worldVertex - camera);
      const float3 origin = camera;
      const float3 end = worldVertex - band * direction;
      const float dist = length(end - origin); 
      float3 step = direction*stepsize[step_phase];

      int3 voxel;
      float3 voxelPos = origin;
      float travelled = 0;
      for(; travelled < dist; travelled += stepsize[step_phase]){
        float3 voxelScaled = floorf(voxelPos * inverseVoxelSize);
        if((voxelScaled.x < size) && (voxelScaled.y < size) &&
            (voxelScaled.z < size) && (voxelScaled.x >= 0) &&
            (voxelScaled.y >= 0) && (voxelScaled.z >= 0)){
          voxel = make_int3(voxelScaled);
          if(!map_index.fetch_octant(voxel.x, voxel.y, voxel.z, tree_depth[step_phase])){
            uint k = map_index.hash(voxel.x, voxel.y, voxel.z);
            unsigned int idx = ++(voxelCount[step_phase]);
            if(idx < reserved[step_phase]) {
              allocationList[step_phase][idx] = k;
            }
          }
        }
        if(step_phase == 0 && 
            dist - travelled < (3 * stepsize[0])) {
          step_phase = 1;
          step = direction*stepsize[step_phase];
        }
        voxelPos +=step;
      }
    }
  }
  written[0] = voxelCount[0];
  written[1] = voxelCount[1];
}
#endif
