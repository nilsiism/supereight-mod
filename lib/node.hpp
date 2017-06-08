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

#ifndef NODE_H
#define NODE_H

#include <commons.h>
#include <cutil_math.h>
#include <atomic>
#include <octree_defines.h>
#include <memory_pool.hpp>

template <typename T>
class Node {

public:

  Node(){
    for (unsigned int i = 0; i < 8; i++)
      child_ptr_[i] = NULL;
    }

    virtual ~Node(){};

    Node *& child(const uint x, const uint y, const uint z);
    Node *& child(int offset );

    static int size(){ return sizeof(Node *) + sizeof(uint3) + 
                                     sizeof(int) + sizeof(short2) + 
                                     sizeof(int); };

    virtual bool isLeaf(){ return false; }

protected:
    Node *child_ptr_[8];
};

template <typename T>
class VoxelBlock: public Node<T> {

  public:

    typedef kfusion_voxel_traits<T> traits_type;
    typedef typename traits_type::ComputeType compute_type;
    typedef typename traits_type::StoredType stored_type;
    static constexpr unsigned int side = BLOCK_SIDE;
    static constexpr unsigned int sideSq = side*side;

    static constexpr compute_type empty() { 
      return traits_type::empty(); 
    }
    static constexpr stored_type initValue() { 
      return traits_type::initValue();
    }
    stored_type translate(const compute_type value) { 
      return traits_type::translate(value); 
    }
    compute_type translate(const stored_type value) const {
      return traits_type::translate(value);
    }


    // Default constructor -- initialises voxel data to their initial values.
    VoxelBlock();

    bool isLeaf(){ return true; }

    uint3 coordinates(){ return coordinates_; }
    void coordinates(const uint3 c){ coordinates_ = c; }

    compute_type data(const uint3 pos) const;
    void data(const uint3 pos, compute_type& value);

    inline float vs3(const uint x, const uint y, const uint z) const {
      return voxel_block_[(x) +
        (y)*side +
        (z)*sideSq].x;
    }

    void last_integrated_frame(const int frame){ 
      last_integrated_frame_ = frame;
    }

    int last_integrated_frame() const {
      return last_integrated_frame_;
    }

    void active(const bool a){ active_ = a; }
    bool active(){ return active_; }

    stored_type * getBlockRawPtr(){ return voxel_block_; }
    static constexpr int size(){ return sizeof(VoxelBlock<T>); }

  private:
    VoxelBlock(const VoxelBlock&) = delete;
    uint3 coordinates_;
    stored_type voxel_block_[side*sideSq]; // Brick of data.
    int last_integrated_frame_;
    bool active_;
};

template <typename T>
VoxelBlock<T>:: VoxelBlock(){
    last_integrated_frame_ = 0;
    coordinates_ = make_uint3(0);
    for (unsigned int i = 0; i < side*sideSq; i++){
      voxel_block_[i] = initValue();
    }
  }

template <typename T>
inline Node<T> *& Node<T>::child(const uint x, const uint y, const uint z){
  return (child_ptr_[x + y*2 + z*4]);
}

template <typename T>
inline Node<T> *& Node<T>::child(int offset){
  return child_ptr_[offset];
}

template <typename T>
inline typename VoxelBlock<T>::compute_type VoxelBlock<T>::data(const uint3 pos) const {
  uint3 offset = pos - coordinates_;
  const stored_type& data = voxel_block_[offset.x + offset.y*side + 
                                         offset.z*sideSq];
  return translate(data);
}

template <typename T>
inline void VoxelBlock<T>::data(const uint3 pos, 
                               compute_type &value){
  uint3 offset = pos - coordinates_;
  voxel_block_[offset.x + offset.y*side + offset.z*sideSq] = translate(value);
}
#endif
