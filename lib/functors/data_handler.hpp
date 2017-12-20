#ifndef DATA_HANDLER_HPP
#define DATA_HANDLER_HPP
#include "math_utils.h"
#include "node.hpp"

template <typename SpecialisedHandlerT, typename NodeT>
class DataHandlerBase {
  typename NodeT::compute_type get() {
    return static_cast<SpecialisedHandlerT *>(this)->get();
  } 
  void set(const typename NodeT::compute_type& val) {
    static_cast<SpecialisedHandlerT *>(this)->set(val);
  }
};

template<typename FieldType>
class VoxelBlockHandler : 
  DataHandlerBase<VoxelBlockHandler<FieldType>, VoxelBlock<FieldType> > {

public:
  VoxelBlockHandler(VoxelBlock<FieldType>* ptr, int3 v) : 
    _block(ptr), _voxel(v) {}

  typename VoxelBlock<FieldType>::compute_type get() {
    return _block->data(_voxel);
  }

  void set(const typename VoxelBlock<FieldType>::compute_type& val) {
    _block->data(_voxel, val);
  }

  private:
    VoxelBlock<FieldType> * _block;  
    int3 _voxel;
};

template<typename FieldType>
class NodeHandler: DataHandlerBase<NodeHandler<FieldType>, Node<FieldType> > {
  public:
    NodeHandler(Node<FieldType>* ptr) : _node(ptr) {}

    typename Node<FieldType>::compute_type get() {
      return _node->value_;
    }

    void set(const typename Node<FieldType>::compute_type& val) {
      _node->value_ = val;
    }

  private:
    Node<FieldType> * _node; 
};


#endif
