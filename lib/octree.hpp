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

#ifndef OCTREE_H
#define OCTREE_H

#include <vector_types.h>
#include <octree_defines.h>
#include <utils/morton_utils.hpp>
#include <cutil_math.h>
#include <commons.h>
#include <algorithm>
#include <parallel/algorithm>
#include <node.hpp>
#include <memory_pool.hpp>
#include <algorithms/mapping.hpp>

#define MAX_BITS 10
#define CAST_STACK_DEPTH 23


uint MASK[] = {
  0x38000000, // 111 000 000 000 000 000 000 000 000 000
  0x3F000000, // 111 111 000 000 000 000 000 000 000 000
  0x3FE00000, // 111 111 111 000 000 000 000 000 000 000
  0x3FFC0000, // 111 111 111 111 000 000 000 000 000 000
  0x3FFF8000, // 111 111 111 111 111 000 000 000 000 000
  0x3FFFF000, // 111 111 111 111 111 111 000 000 000 000
  0x3FFFFE00, // 111 111 111 111 111 111 111 000 000 000
  0x3FFFFFC0, // 111 111 111 111 111 111 111 111 000 000
  0x3FFFFFF8, // 111 111 111 111 111 111 111 111 111 000
  0x3FFFFFFF  // 111 111 111 111 111 111 111 111 111 111
};

inline int __float_as_int(float value){

  union float_as_int {
    float f;
    int i;
  };

  float_as_int u;
  u.f = value;
  return u.i;
}

inline float __int_as_float(int value){

  union int_as_float {
    int i;
    float f;
  };

  int_as_float u;
  u.i = value;
  return u.f;
}

template <typename T>
class Octree
{

public:

  typedef kfusion_voxel_traits<T> traits_type;
  typedef typename traits_type::ComputeType compute_type;
  typedef typename traits_type::StoredType stored_type;
  compute_type empty() const { return traits_type::empty(); }

  Octree(){
  };

  ~Octree(){
    deallocateTree();
  }

  void init(int size, float3 dim);

  inline int size(){ return size_; }
  inline float3 dim() const { return dim_; }
  float3 pos(const uint3 p) const;

  unsigned int hash(const int x, const int y, const int z) {
    constexpr int morton_mask = MAX_BITS - log2(blockSide) - 1;
    return compute_morton(make_uint3(x, y, z)) & MASK[morton_mask];   
  }
  
  compute_type get(const int x, const int y, const int z) const;
  compute_type get(const float3 pos) const;
  Aggregate<T> * fetch(const int x, const int y, const int z) const;
  Aggregate<T> * fetch(const float x, const float y, const float z) const;

  float interp(const float3 pos) const;
  float3 grad(const float3 pos) const;


  float4 raycast(const uint2 pos, const Matrix4 view,
      const float nearPlane, const float farPlane, const float mu,
      const float step, const float largeStep, int frame);

  // Returns a list of allocated blocks. The active switch indicates whether
  // the blocks should be visible or not.
  void getBlockList(std::vector<Aggregate<T> *>& blocklist, bool active);

  // Allocate the list of voxels, passed as a list of morton numbers.
  bool allocate(uint *keys, int num_elem);

  int leavesCount();
  int nodeCount();

  void integrateFrame(const Matrix4 &pose, const float4& k,
      const float *depthmap, const uint2 &imageSize,
      const float mu, const int frame);

  void printMemStats(){
    // memory.printStats();
  };
  bool print_kernel_timing;

private:

  struct stack_entry {
    int scale;
    Node<T> * parent;
    float t_max;
  };

  Node<T> * root_;
  int size_;
  float3 dim_;
  int max_level_;
  MemoryPool<Aggregate<T> > block_memory_;

  // Compile-time constant expressions
  // # of voxels per side in a voxel block
  static constexpr unsigned int blockSide = BLOCK_SIDE;
  // maximum tree depth in bits
  static constexpr unsigned int max_depth = ((sizeof(morton_type)*8)/3);
  // Tree depth at which blocks are found
  static constexpr unsigned int block_depth = max_depth -
    static_cast<unsigned int>(log2(BLOCK_SIDE));

  // Allocation specific variables
  uint * allocationList_;
  uint * keys_at_level_;
  int reserved_;

  // Camera parameters
  float maxweight_;  // maximum weight

  // Private implementation of cached methods
  compute_type get(const int x, const int y, const int z, Aggregate<T>* cached) const;
  compute_type get(const float3 pos, Aggregate<T>* cached) const;
  float vs(const uint x, const uint y, const uint z, Aggregate<T>*& start) const;

  // Build a list of voxels to be allocated given the current input frame.
  // Returns the number of voxels to be allocated.
  unsigned int buildAllocationList(const Matrix4 &pose, const float4& k,
      const float *depthmap, const uint2 &imageSize,
      const float mu, const int frame);


  // Parallel allocation of a given tree level for a set of input keys.
  // Pre: levels above target_level must have been already allocated
  bool allocateLevel(uint * keys, int num_tasks, int target_level);

  // Masks code with the appropriate bitmask for the input three level.
  unsigned int getMortonAtLevel(uint code, int level);

  void reserveBuffers(const int n);
  bool getKeysAtLevel(const uint * inputKeys, uint *outpuKeys, int num_keys, int level);
  uint3 getChildFromCode(int code, int level);

  // General helpers

  int leavesCountRecursive(Node<T> *);
  int nodeCountRecursive(Node<T> *);
  void getActiveBlockList(Node<T> *, std::vector<Aggregate<T> *>& blocklist);
  void getAllocatedBlockList(Node<T> *, std::vector<Aggregate<T> *>& blocklist);

  void deleteNode(Node<T> ** node);
  void deallocateTree(){ deleteNode(&root_); }
};

template <typename T>
inline typename Octree<T>::compute_type Octree<T>::get(const float3 p) const {

  const uint3 pos = make_uint3(p.x, p.y, p.z);

  Node<T> * n = root_;
  if(!n) {
    return empty();
  }

  // Get the block.

  uint edge = size_ >> 1;
  for(; edge >= blockSide; edge = edge >> 1){
    n = n->child((pos.x & edge) > 0u, (pos.y & edge) > 0u, (pos.z & edge) > 0u);
    if(!n){
    return empty();
    }
  }

  // Get the element in the voxel block
  return static_cast<Aggregate<T>*>(n)->data(pos);
}

template <typename T>
inline typename Octree<T>::compute_type Octree<T>::get(const float3 p, 
    Aggregate<T>* cached) const {

  const uint3 pos = make_uint3((p.x * size_ / dim_.x),
                               (p.y * size_ / dim_.y),
                               (p.z * size_ / dim_.z));

  if(cached != NULL){
    uint3 lower = cached->coordinates();
    uint3 upper = lower + (blockSide-1);
    if(in(pos.x, lower.x, upper.x) && in(pos.y, lower.y, upper.y) &&
       in(pos.z, lower.z, upper.z)){
      return cached->data(pos);
    }
  }

  Node<T> * n = root_;
  if(!n) {
    return empty();
  }

  // Get the block.

  uint edge = size_ >> 1;
  for(; edge >= blockSide; edge = edge >> 1){
    n = n->child((pos.x & edge) > 0u, (pos.y & edge) > 0u, (pos.z & edge) > 0u);
    if(!n){
    return empty();
    }
  }

  // Get the element in the voxel block
  return static_cast<Aggregate<T>*>(n)->data(pos);
}

template <typename T>
inline typename Octree<T>::compute_type Octree<T>::get(const int x,
    const int y, const int z) const {
  Node<T> * n = root_;
  int level = 1;

  if(!n) {
    return empty();
  }

  const uint ux = x;
  const uint uy = y;
  const uint uz = z;
  uint edge = size_ >> 1;
  for(; edge >= blockSide; edge = edge >> 1){
    n = n->child((ux & edge) > 0u, (uy & edge) > 0u, (uz & edge) > 0u);
    if(!n){
      return empty();
    }
  }

  return static_cast<Aggregate<T> *>(n)->data(make_uint3(x, y, z));
}

template <typename T>
inline typename Octree<T>::compute_type Octree<T>::get(const int x,
   const int y, const int z, Aggregate<T>* cached) const {

  const uint ux = x;
  const uint uy = y;
  const uint uz = z;

  if(cached != NULL){
    const uint3 lower = cached->coordinates();
    const uint3 upper = lower + (blockSide-1);
    if(in(x, lower.x, upper.x) && in(y, lower.y, upper.y) &&
       in(z, lower.z, upper.z)){
      return cached->data(make_uint3(x, y, z));
    }
  }

  Node<T> * n = root_;
  int level = 1;

  if(!n) {
    return empty();
  }

  uint edge = size_ >> 1;
  for(; edge >= blockSide; edge = edge >> 1){
    n = n->child((ux & edge) > 0u, (uy & edge) > 0u, (uz & edge) > 0u);
    if(!n){
      return empty();
    }
  }

  return static_cast<Aggregate<T> *>(n)->data(make_uint3(ux, uy, uz));
}

template <typename T>
void Octree<T>::deleteNode(Node<T> **node){

  if(*node){
    for (int i = 0; i < 8; i++) {
      if((*node)->child(i)){
        deleteNode(&(*node)->child(i));
      }
    }
    if(!(*node)->isLeaf()){
      delete *node;
      *node = NULL;
    }
  }
}


template <typename T>
void Octree<T>::init(int size, float3 dim) {
  size_ = size;
  dim_ = dim;
  max_level_ = log2(size);
  root_ = new Node<T>();
  // root_->edge(size_);
  reserved_ = 0;
  maxweight_ = maxweight; // from constant_parameters.h
}

template <typename T>
float3 Octree<T>::pos(const uint3 p) const {
  return make_float3((p.x + 0.5f) * dim_.x / size_,
      (p.y + 0.5f) * dim_.y / size_, (p.z + 0.5f) * dim_.z / size_);
}

template <typename T>
inline Aggregate<T> * Octree<T>::fetch(const float x, const float y, 
   const float z) const {
  return fetch((int)x, (int)y, (int)z);
}

template <typename T>
inline Aggregate<T> * Octree<T>::fetch(const int x, const int y, 
   const int z) const {

  Node<T> * n = root_;
  if(!n) {
    return NULL;
  }

  // Get the block.
  uint edge = size_ / 2;
  for(; edge >= blockSide; edge /= 2){
    n = n->child((x & edge) > 0u, (y & edge) > 0u, (z & edge) > 0u);
    if(!n){
      return NULL;
    }
  }
  return static_cast<Aggregate<T>* > (n);
}

inline bool operator<=(const uint3 a, const int3 b){
  return((a.x <= b.x) && (a.y <= b.y) && (a.z <= b.z));
}

template <typename T>
float Octree<T>::interp(const float3 pos) const {

  const int3 base = make_int3(floorf(pos));
  const float3 factor = fracf(pos);
  const int3 lower = max(base, make_int3(0));
  const int3 upper = min(base + make_int3(1),
      make_int3(size_) - make_int3(1));

  Aggregate<T> * n = fetch(lower.x, lower.y, lower.z);
  if(n){
    const uint3 ul = n->coordinates() + Aggregate<T>::side;
    // Local interpolation
    if( upper.x < ul.x && upper.y < ul.y && upper.z < ul.z){

      stored_type* data = n->getBlockRawPtr();
      const uint3 offset = n->coordinates();
      uint3 lower_offset = make_uint3(lower.x - offset.x,
          lower.y - offset.y,
          lower.z - offset.z);

      uint3 upper_offset = make_uint3(upper.x - offset.x,
          upper.y - offset.y,
          upper.z - offset.z);

      static constexpr uint edge = blockSide;
      static constexpr uint edge2 = edge * edge;

      float value = (((data[lower_offset.x + lower_offset.y*edge + lower_offset.z*edge2].x * (1 - factor.x)
              + data[upper_offset.x + lower_offset.y*edge + lower_offset.z*edge2].x * factor.x) * (1 - factor.y)
            + (data[lower_offset.x + (upper_offset.y)*edge + lower_offset.z*edge2].x * (1 - factor.x)
              + data[upper_offset.x + (upper_offset.y)*edge + lower_offset.z*edge2].x * factor.x) * factor.y)
          * (1 - factor.z)
          + ((data[lower_offset.x + lower_offset.y*edge + (upper_offset.z)*edge2].x * (1 - factor.x)
              + data[upper_offset.x + lower_offset.y*edge + (upper_offset.z)*edge2].x * factor.x)
            * (1 - factor.y)
            + (data[lower_offset.x + (upper_offset.y)*edge + (upper_offset.z)*edge2].x * (1 - factor.x)
              + data[upper_offset.x + (upper_offset.y)*edge + (upper_offset.z)*edge2].x * factor.x)
            * factor.y) * factor.z) * 0.00003051944088f;
      return value;
    }
  }

  return (((get(lower.x, lower.y, lower.z, n).x * (1 - factor.x)
          + get(upper.x, lower.y, lower.z, n).x * factor.x) * (1 - factor.y)
          + ( get(lower.x, upper.y, lower.z, n).x * (1 - factor.x)
          + get(upper.x, upper.y, lower.z, n).x * factor.x) * factor.y)
          * (1 - factor.z)
          + ((  get(lower.x, lower.y, upper.z, n).x * (1 - factor.x)
          + get(upper.x, lower.y, upper.z, n).x * factor.x)
          * (1 - factor.y)
          + ( get(lower.x, upper.y, upper.z, n).x * (1 - factor.x)
          + get(upper.x, upper.y, upper.z, n).x * factor.x)
          * factor.y) * factor.z);
}


template <typename T>
float3 Octree<T>::grad(const float3 pos) const {

   int3 base = make_int3(floorf(pos));
   float3 factor = fracf(pos);
   int3 lower_lower = max(base - make_int3(1), make_int3(0));
   int3 lower_upper = max(base, make_int3(0));
   int3 upper_lower = min(base + make_int3(1),
      make_int3(size_) - make_int3(1));
   int3 upper_upper = min(base + make_int3(2),
      make_int3(size_) - make_int3(1));
   int3 & lower = lower_upper;
   int3 & upper = upper_lower;

  const int lower_code = compute_morton(lower_lower) & MASK[block_depth + MAX_BITS - max_level_ - 1];
  const int upper_code = compute_morton(upper_upper) & MASK[block_depth + MAX_BITS - max_level_ - 1];

  float3 gradient;

  if (!(lower_code ^ upper_code)){
//  std::cout << "Doing local gradient here!" << std::endl;
    Aggregate<T> * n = fetch(lower_lower.x, lower_lower.y, lower_lower.z);
    if(n != NULL){
      int3 offset = make_int3(n->coordinates());
      lower_lower = lower_lower - offset;
      lower_upper = lower_upper - offset;
      upper_lower = upper_lower - offset;
      upper_upper = upper_upper - offset;

      gradient.x = (((n->vs3(upper_lower.x, lower.y, lower.z)
                    - n->vs3(lower_lower.x, lower.y, lower.z)) * (1 - factor.x)
                   + (n->vs3(upper_upper.x, lower.y, lower.z)
                    - n->vs3(lower_upper.x, lower.y, lower.z)) * factor.x)
                    * (1 - factor.y)
                    + ((n->vs3(upper_lower.x, upper.y, lower.z)
                        - n->vs3(lower_lower.x, upper.y, lower.z)) * (1 - factor.x)
                      + (n->vs3(upper_upper.x, upper.y, lower.z)
                        - n->vs3(lower_upper.x, upper.y, lower.z))
                      * factor.x) * factor.y) * (1 - factor.z)
        + (((n->vs3(upper_lower.x, lower.y, upper.z)
                - n->vs3(lower_lower.x, lower.y, upper.z)) * (1 - factor.x)
              + (n->vs3(upper_upper.x, lower.y, upper.z)
                - n->vs3(lower_upper.x, lower.y, upper.z))
              * factor.x) * (1 - factor.y)
            + ((n->vs3(upper_lower.x, upper.y, upper.z)
                - n->vs3(lower_lower.x, upper.y, upper.z))
              * (1 - factor.x)
              + (n->vs3(upper_upper.x, upper.y, upper.z)
                - n->vs3(lower_upper.x, upper.y, upper.z))
              * factor.x) * factor.y) * factor.z;

  gradient.y = (((n->vs3(lower.x, upper_lower.y, lower.z)
          - n->vs3(lower.x, lower_lower.y, lower.z)) * (1 - factor.x)
        + (n->vs3(upper.x, upper_lower.y, lower.z)
          - n->vs3(upper.x, lower_lower.y, lower.z)) * factor.x)
      * (1 - factor.y)
      + ((n->vs3(lower.x, upper_upper.y, lower.z)
          - n->vs3(lower.x, lower_upper.y, lower.z)) * (1 - factor.x)
        + (n->vs3(upper.x, upper_upper.y, lower.z)
          - n->vs3(upper.x, lower_upper.y, lower.z))
        * factor.x) * factor.y) * (1 - factor.z)
    + (((n->vs3(lower.x, upper_lower.y, upper.z)
            - n->vs3(lower.x, lower_lower.y, upper.z)) * (1 - factor.x)
          + (n->vs3(upper.x, upper_lower.y, upper.z)
            - n->vs3(upper.x, lower_lower.y, upper.z))
          * factor.x) * (1 - factor.y)
        + ((n->vs3(lower.x, upper_upper.y, upper.z)
            - n->vs3(lower.x, lower_upper.y, upper.z))
          * (1 - factor.x)
          + (n->vs3(upper.x, upper_upper.y, upper.z)
            - n->vs3(upper.x, lower_upper.y, upper.z))
          * factor.x) * factor.y) * factor.z;

  gradient.z = (((n->vs3(lower.x, lower.y, upper_lower.z)
          - n->vs3(lower.x, lower.y, lower_lower.z)) * (1 - factor.x)
        + (n->vs3(upper.x, lower.y, upper_lower.z)
          - n->vs3(upper.x, lower.y, lower_lower.z)) * factor.x)
      * (1 - factor.y)
      + ((n->vs3(lower.x, upper.y, upper_lower.z)
          - n->vs3(lower.x, upper.y, lower_lower.z)) * (1 - factor.x)
        + (n->vs3(upper.x, upper.y, upper_lower.z)
          - n->vs3(upper.x, upper.y, lower_lower.z))
        * factor.x) * factor.y) * (1 - factor.z)
    + (((n->vs3(lower.x, lower.y, upper_upper.z)
            - n->vs3(lower.x, lower.y, lower_upper.z)) * (1 - factor.x)
          + (n->vs3(upper.x, lower.y, upper_upper.z)
            - n->vs3(upper.x, lower.y, lower_upper.z))
          * factor.x) * (1 - factor.y)
        + ((n->vs3(lower.x, upper.y, upper_upper.z)
            - n->vs3(lower.x, upper.y, lower_upper.z))
          * (1 - factor.x)
          + (n->vs3(upper.x, upper.y, upper_upper.z)
            - n->vs3(upper.x, upper.y, lower_upper.z))
          * factor.x) * factor.y) * factor.z;

  return gradient
    * make_float3(dim_.x / size_, dim_.y / size_, dim_.z / size_)
    * (0.5f * 0.00003051944088f);
    }
  }

  Aggregate<T> * n = fetch(base.x, base.y, base.z);

  gradient.x = (((get(upper_lower.x, lower.y, lower.z, n).x
          - get(lower_lower.x, lower.y, lower.z, n).x) * (1 - factor.x)
        + (get(upper_upper.x, lower.y, lower.z, n).x
          - get(lower_upper.x, lower.y, lower.z, n).x) * factor.x)
      * (1 - factor.y)
      + ((get(upper_lower.x, upper.y, lower.z, n).x
          - get(lower_lower.x, upper.y, lower.z, n).x) * (1 - factor.x)
        + (get(upper_upper.x, upper.y, lower.z, n).x
          - get(lower_upper.x, upper.y, lower.z, n).x)
        * factor.x) * factor.y) * (1 - factor.z)
    + (((get(upper_lower.x, lower.y, upper.z, n).x
            - get(lower_lower.x, lower.y, upper.z, n).x) * (1 - factor.x)
          + (get(upper_upper.x, lower.y, upper.z, n).x
            - get(lower_upper.x, lower.y, upper.z, n).x)
          * factor.x) * (1 - factor.y)
        + ((get(upper_lower.x, upper.y, upper.z, n).x
            - get(lower_lower.x, upper.y, upper.z, n).x)
          * (1 - factor.x)
          + (get(upper_upper.x, upper.y, upper.z, n).x
            - get(lower_upper.x, upper.y, upper.z, n).x)
          * factor.x) * factor.y) * factor.z;

  gradient.y = (((get(lower.x, upper_lower.y, lower.z, n).x
          - get(lower.x, lower_lower.y, lower.z, n).x) * (1 - factor.x)
        + (get(upper.x, upper_lower.y, lower.z, n).x
          - get(upper.x, lower_lower.y, lower.z, n).x) * factor.x)
      * (1 - factor.y)
      + ((get(lower.x, upper_upper.y, lower.z, n).x
          - get(lower.x, lower_upper.y, lower.z, n).x) * (1 - factor.x)
        + (get(upper.x, upper_upper.y, lower.z, n).x
          - get(upper.x, lower_upper.y, lower.z, n).x)
        * factor.x) * factor.y) * (1 - factor.z)
    + (((get(lower.x, upper_lower.y, upper.z, n).x
            - get(lower.x, lower_lower.y, upper.z, n).x) * (1 - factor.x)
          + (get(upper.x, upper_lower.y, upper.z, n).x
            - get(upper.x, lower_lower.y, upper.z, n).x)
          * factor.x) * (1 - factor.y)
        + ((get(lower.x, upper_upper.y, upper.z, n).x
            - get(lower.x, lower_upper.y, upper.z, n).x)
          * (1 - factor.x)
          + (get(upper.x, upper_upper.y, upper.z, n).x
            - get(upper.x, lower_upper.y, upper.z, n).x)
          * factor.x) * factor.y) * factor.z;

  gradient.z = (((get(lower.x, lower.y, upper_lower.z, n).x
          - get(lower.x, lower.y, lower_lower.z, n).x) * (1 - factor.x)
        + (get(upper.x, lower.y, upper_lower.z, n).x
          - get(upper.x, lower.y, lower_lower.z, n).x) * factor.x)
      * (1 - factor.y)
      + ((get(lower.x, upper.y, upper_lower.z, n).x
          - get(lower.x, upper.y, lower_lower.z, n).x) * (1 - factor.x)
        + (get(upper.x, upper.y, upper_lower.z, n).x
          - get(upper.x, upper.y, lower_lower.z, n).x)
        * factor.x) * factor.y) * (1 - factor.z)
    + (((get(lower.x, lower.y, upper_upper.z, n).x
            - get(lower.x, lower.y, lower_upper.z, n).x) * (1 - factor.x)
          + (get(upper.x, lower.y, upper_upper.z, n).x
            - get(upper.x, lower.y, lower_upper.z, n).x)
          * factor.x) * (1 - factor.y)
        + ((get(lower.x, upper.y, upper_upper.z, n).x
            - get(lower.x, upper.y, lower_upper.z, n).x)
          * (1 - factor.x)
          + (get(upper.x, upper.y, upper_upper.z, n).x
            - get(upper.x, upper.y, lower_upper.z, n).x)
          * factor.x) * factor.y) * factor.z;

  return gradient
    * make_float3(dim_.x / size_, dim_.y / size_, dim_.z / size_)
    * (0.5f);

}

template <typename T>
int Octree<T>::leavesCount(){
  return leavesCountRecursive(root_);
}

template <typename T>
int Octree<T>::leavesCountRecursive(Node<T> * n){

  if(!n) return 0;

  if(n->isLeaf()){
    return 1;
  }

  int sum = 0;

  for (int i = 0; i < 8; i++){
    sum += leavesCountRecursive(n->child(i));
  }

  return sum;
}

template <typename T>
int Octree<T>::nodeCount(){
  return nodeCountRecursive(root_);
}

template <typename T>
int Octree<T>::nodeCountRecursive(Node<T> * node){
  if (!node) {
    return 0;
  }

  int n = 1;
  for (int i = 0; i < 8; ++i) {
    n += (n ? nodeCountRecursive((node)->child(i)) : 0);
  }
  return n;
}

template <typename T>
inline uint Octree<T>::getMortonAtLevel(uint code, int level){
  const int shift = MAX_BITS - max_level_;
  return code & MASK[level + shift-1];
}

template <typename T>
inline bool Octree<T>::getKeysAtLevel(const uint * inputKeys, uint * outputKeys,  int num_keys, int level){

  const int shift = MAX_BITS - max_level_;
  outputKeys[0] = inputKeys[0] & MASK[level + shift-1];
  for (unsigned int i = 1; i < num_keys; i++){
    outputKeys[i] = inputKeys[i] & MASK[level + shift-1];
  }

  return true;
}

template <typename T>
void Octree<T>::reserveBuffers(const int n){

  if(n > reserved_){
    // std::cout << "Reserving " << n << " entries in allocation buffers" << std::endl;
    delete[] keys_at_level_;
    keys_at_level_ = new uint[n];
    reserved_ = n;
  }
  block_memory_.reserve(n);
}

template <typename T>
bool Octree<T>::allocate(uint *keys, int num_elem){

#ifdef _OPENMP
  __gnu_parallel::sort(keys, keys+num_elem);
#else
std::sort(keys, keys+num_elem);
#endif

  num_elem = unique(keys, num_elem);
  reserveBuffers(num_elem);

  int last_elem = 0;
  const int leaf_level = max_level_ - log2(blockSide);
  for (int level = 1; level <= leaf_level; level++){
    getKeysAtLevel(keys, keys_at_level_, num_elem, level);
    last_elem = unique(keys_at_level_, num_elem);
    bool success = allocateLevel(keys_at_level_, last_elem, level);
  }
  return true;
}

template <typename T>
bool Octree<T>::allocateLevel(uint * keys, int num_tasks, int target_level){

  int leaves_level = max_level_ - log2(blockSide);

#pragma omp parallel for
  for (int i = 0; i < num_tasks; i++){
    Node<T> ** n = &root_;
    Node<T> * parent;
    int myKey = keys[i];
    int edge = size_/2;

    for (int level = 1; level <= target_level; ++level){

      uint3 child = getChildFromCode(myKey, level);
      int index = child.x + child.y*2 + child.z*4;
      parent = *n;
      n = &(*n)->child(index);

      if(!(*n)){
        if(level == leaves_level){
          *n = block_memory_.acquire_block();
          static_cast<Aggregate<T> *>(*n)->coordinates(unpack_morton(myKey));
          static_cast<Aggregate<T> *>(*n)->active(true);
        }
        else  *n = new Node<T>();
      }
      edge /= 2;
    }
  }
  return true;
}



template <typename T>
inline uint3 Octree<T>::getChildFromCode(int code, int level){

  int shift = max_level_ - level;
  code = code >> shift*3;

  uint3 coordinates = make_uint3(code & 0x01, (code >> 1) & 0x01, (code >> 2) & 0x01);

  return coordinates;

}

/**
 * A modified version of the ray-caster introduced in the paper:
 * https://research.nvidia.com/publication/efficient-sparse-voxel-octrees
 *
 * Original code available at:
 * https://code.google.com/p/efficient-sparse-voxel-octrees/
 *
**/

template <typename T>
float4 Octree<T>::raycast(const uint2 position, const Matrix4 view,
    const float nearPlane, const float farPlane, const float mu,
    const float step, const float largeStep, int frame){

  const float3 origin = get_translation(view);
  float3 direction = rotate(view, make_float3(position.x, position.y, 1.f));
  struct stack_entry stack[CAST_STACK_DEPTH];
  static const float epsilon = exp2f(-log2(size_));

  if(fabsf(direction.x) < epsilon) direction.x = copysignf(epsilon, direction.x);
  if(fabsf(direction.y) < epsilon) direction.y = copysignf(epsilon, direction.y);
  if(fabsf(direction.z) < epsilon) direction.z = copysignf(epsilon, direction.z);

  // Scaling the origin to resides between coordinates [1,2]
  const float3 scaled_origin = origin/dim_ + 1;

  // Precomputing the coefficients of tx(x), ty(y) and tz(z)
  // The octree is assumed to reside at coordinates [1,2]

  float3 t_coef = -1/fabs(direction);
  float3 t_bias = t_coef * scaled_origin;

  // Build the octanct mask to mirror the coordinate system
  // so that each ray direction component is negative.

  int octant_mask = 7;
  if(direction.x > 0.0f) octant_mask ^=1, t_bias.x = 3.0f * t_coef.x - t_bias.x;
  if(direction.y > 0.0f) octant_mask ^=2, t_bias.y = 3.0f * t_coef.y - t_bias.y;
  if(direction.z > 0.0f) octant_mask ^=4, t_bias.z = 3.0f * t_coef.z - t_bias.z;

  // Find the active t-span

  float t_min = fmaxf(fmaxf(2.0f * t_coef.x - t_bias.x, 2.0f * t_coef.y - t_bias.y), 2.0f * t_coef.z - t_bias.z);
  float t_max = fminf(fminf(t_coef.x - t_bias.x, t_coef.y - t_bias.y), t_coef.z - t_bias.z);
  float h = t_max;
  t_min = fmaxf(t_min, 0.0f);
  t_max = fminf(t_max, farPlane/dim_.x);

  Node<T> * parent = root_;
  Node<T> * child;
  int idx = 0;
  float3 pos = make_float3(1.0f, 1.0f, 1.0f);
  int scale = CAST_STACK_DEPTH-1;
  float scale_exp2 = 0.5f;
  int min_scale = CAST_STACK_DEPTH - log2(size_/blockSide);
  float last_sample = 1.0f;
  float last_t = 0;
  float stepsize = largeStep;// largeStep = 0.75*mu

  if (1.5f * t_coef.x- t_bias.x > t_min) idx ^= 1, pos.x = 1.5f;
  if (1.5f * t_coef.y- t_bias.y > t_min) idx ^= 2, pos.y = 1.5f;
  if (1.5f * t_coef.z- t_bias.z > t_min) idx ^= 4, pos.z = 1.5f;

  while (scale < CAST_STACK_DEPTH) {
    float3 t_corner = pos * t_coef - t_bias;
    float tc_max = fminf(fminf(t_corner.x, t_corner.y), t_corner.z);

    int child_idx = idx ^ octant_mask ^ 7;
    child = parent->child(child_idx);

    if (scale == min_scale && child != NULL){

        // check against near and far plane
        float tnear = t_min  * dim_.x;
        float tfar = tc_max * dim_.x;
        // const float tfar = fminf(fminf(t_corner.x, t_corner.y), t_corner.z);
        if (tnear < tfar) {
          // first walk with largesteps until we found a hit
          float t = last_t;
          float f_t = last_sample;
          float f_tt = last_sample;
          if (f_t > 0.f) {
            for (t = tnear; t < tfar; t += stepsize) {
              f_tt = get(origin + direction * t, static_cast<Aggregate<T>*>(child)).x;
              if(f_tt < 0.1f){
                f_tt = interp(origin + direction * t);
              }
              if (f_tt < 0.f)                  // got it, jump out of inner loop
                break;
              stepsize = fmaxf(f_tt*mu, step);
              f_t = f_tt;
            }
            if (f_tt < 0.f) {           // got it, calculate accurate intersection
              t = t + stepsize * f_tt / (f_t - f_tt);
              float4 hit =  make_float4(origin + direction * t, t);
              return hit;
            }
            last_sample = f_tt;
          }
        }
    } else if (child != NULL && t_min <= t_max){  // If the child is valid, descend the tree hierarchy.

      float tv_max = fminf(t_max, tc_max);
      float half = scale_exp2 * 0.5f;
      float3 t_center = half * t_coef + t_corner;

      // Descend to the first child if the resulting t-span is non-empty.

      if (tc_max < h) {
        stack[scale] = {scale, parent, t_max};
      }

      h = tc_max;
      parent = child;

      idx = 0;
      scale--;
      scale_exp2 = half;
      if (t_center.x > t_min) idx ^= 1, pos.x += scale_exp2;
      if (t_center.y > t_min) idx ^= 2, pos.y += scale_exp2;
      if (t_center.z > t_min) idx ^= 4, pos.z += scale_exp2;

      t_max = tv_max;
      child = NULL;
      continue;
    }

    // ADVANCE
    // Step along the ray.

    int step_mask = 0;

    if (t_corner.x <= tc_max) step_mask ^= 1, pos.x -= scale_exp2;
    if (t_corner.y <= tc_max) step_mask ^= 2, pos.y -= scale_exp2;
    if (t_corner.z <= tc_max) step_mask ^= 4, pos.z -= scale_exp2;

    t_min = tc_max;
    idx ^= step_mask;

    // POP if bits flips disagree with ray direction

    if ((idx & step_mask) != 0) {

      // Get the different bits for each component.

      unsigned int differing_bits = 0;
      if ((step_mask & 1) != 0) differing_bits |= __float_as_int(pos.x) ^ __float_as_int(pos.x + scale_exp2);
      if ((step_mask & 2) != 0) differing_bits |= __float_as_int(pos.y) ^ __float_as_int(pos.y + scale_exp2);
      if ((step_mask & 4) != 0) differing_bits |= __float_as_int(pos.z) ^ __float_as_int(pos.z + scale_exp2);

      // Get the scale at which the two differs. Here's there are different subtlelties related to how fp are stored.
      scale = (__float_as_int((float)differing_bits) >> 23) - 127; // position of the highest bit
      scale_exp2 = __int_as_float((scale - CAST_STACK_DEPTH + 127) << 23); // exp2f(scale - s_max)
      struct stack_entry&  e = stack[scale];
      parent = e.parent;
      t_max = e.t_max;

      // Round cube position and extract child slot index.

      int shx = __float_as_int(pos.x) >> scale;
      int shy = __float_as_int(pos.y) >> scale;
      int shz = __float_as_int(pos.z) >> scale;
      pos.x = __int_as_float(shx << scale);
      pos.y = __int_as_float(shy << scale);
      pos.z = __int_as_float(shz << scale);
      idx  = (shx & 1) | ((shy & 1) << 1) | ((shz & 1) << 2);

      h = 0.0f;
      child = NULL;

    }
  }
  //if( doTock ) TOCK("raycastOctreeIteration", frame);
  return make_float4(0);
}

// Simple ray-caster in metric space

template <typename T>
void Octree<T>::integrateFrame(const Matrix4 &pose, const float4& k,
                            const float *depthmap,
                            const uint2 &imageSize, const float mu,
                            const int frame) {

  maxweight_ = 100;

  std::vector<Aggregate<T> *> active_list;
  getBlockList(active_list, true);

  Aggregate<T> ** list = active_list.data();
  unsigned int num_active = active_list.size();
  integratePass(list, num_active, depthmap, imageSize, dim_.x/size_,
      inverse(pose), getCameraMatrix(k),  mu, maxweight_, frame);
}

template <typename T>
unsigned int Octree<T>::buildAllocationList(const Matrix4 &pose, const float4& k,
    const float *depthmap, const uint2 &imageSize,
    const float mu, const int frame) {

  const float voxelSize = dim_.x/size_;
  const float inverseVoxelSize = 1/voxelSize;
  constexpr int morton_mask = MAX_BITS - log2(blockSide) - 1;

  const Matrix4 invPose = inverse(pose);
  Matrix4 K = getCameraMatrix(k);
  Matrix4 invK = getInverseCameraMatrix(k);
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

      uint3 voxel;
      float3 voxelPos = origin;
      for(int i = 0; i < numSteps; i++){
        float3 voxelScaled = floorf(voxelPos * inverseVoxelSize);
        if((voxelScaled.x < size()) && (voxelScaled.y < size()) &&
           (voxelScaled.z < size()) && (voxelScaled.x >= 0) &&
           (voxelScaled.y >= 0) && (voxelScaled.z >= 0)){
          voxel = make_uint3(voxelScaled.x, voxelScaled.y, voxelScaled.z);
          Aggregate<T> * n = fetch((int) voxel.x, (int) voxel.y, (int) voxel.z);
          if(!n){
            uint k = compute_morton(voxel) & MASK[morton_mask];
            unsigned int idx = ++voxelCount;
            allocationList_[idx] = k;
          }
          else if(n->last_integrated_frame() < frame){
           integrate(n, depthmap, imageSize, voxelSize,
               invPose, K, mu, maxweight_);
           n->last_integrated_frame(frame);
          }
        }
        voxelPos +=step;
      }
    }
  }
  return voxelCount;
}

template <typename T>
void Octree<T>::getBlockList(std::vector<Aggregate<T>*>& blocklist, bool active){
  Node<T> * n = root_;
  if(!n) return;
  if(active) getActiveBlockList(n, blocklist);
  else getAllocatedBlockList(n, blocklist);
}

template <typename T>
void Octree<T>::getActiveBlockList(Node<T> *n,
    std::vector<Aggregate<T>*>& blocklist){
  using tNode = Node<T>;
  if(!n) return;
  std::queue<tNode *> q;
  q.push(n);
  while(!q.empty()){
    tNode* node = q.front();
    q.pop();

    if(node->isLeaf()){
      Aggregate<T>* block = static_cast<Aggregate<T> *>(node);
      if(block->active()) blocklist.push_back(block);
      continue;
    }

    for(int i = 0; i < 8; ++i){
      if(node->child(i)) q.push(node->child(i));
    }
  }
}

template <typename T>
void Octree<T>::getAllocatedBlockList(Node<T> *n,
    std::vector<Aggregate<T>*>& blocklist){
  using tNode = Node<T>;
  if(!n) return;
  std::queue<tNode *> q;
  q.push(n);
  while(!q.empty()){
    tNode* node = q.front();
    q.pop();

    if(node->isLeaf()){
      Aggregate<T>* block = static_cast<Aggregate<T> *>(n);
       blocklist.push_back(block);
      continue;
    }

    for(int i = 0; i < 8; ++i){
      if(node->child(i)) q.push(node->child(i));
    }
  }

}

#endif // OCTREE_H
