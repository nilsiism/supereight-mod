#ifndef DATA_HANDLER_HPP
#define DATA_HANDLER_HPP
#include "../utils/se_common.h"
#include "../node.hpp"

template <typename SpecialisedHandlerT, typename NodeT>
class DataHandlerBase {
  typename NodeT::value_type get() {
    return static_cast<SpecialisedHandlerT *>(this)->get();
  } 
  void set(const typename NodeT::value_type& val) {
    static_cast<SpecialisedHandlerT *>(this)->set(val);
  }
};

template<typename FieldType>
class VoxelBlockHandler : 
  DataHandlerBase<VoxelBlockHandler<FieldType>, se::VoxelBlock<FieldType> > {

public:
  VoxelBlockHandler(se::VoxelBlock<FieldType>* ptr, Eigen::Vector3i v) : 
    _block(ptr), _voxel(v) {}

  typename se::VoxelBlock<FieldType>::value_type get() {
    return _block->data(_voxel);
  }

  void set(const typename se::VoxelBlock<FieldType>::value_type& val) {
    _block->data(_voxel, val);
  }

  private:
    se::VoxelBlock<FieldType> * _block;  
    Eigen::Vector3i _voxel;
};

template<typename FieldType>
class NodeHandler: DataHandlerBase<NodeHandler<FieldType>, se::Node<FieldType> > {
  public:
    NodeHandler(se::Node<FieldType>* ptr, int i) : _node(ptr), _idx(i) {}

    typename se::Node<FieldType>::value_type get() {
      return _node->value_[_idx];
    }

    void set(const typename se::Node<FieldType>::value_type& val) {
      _node->value_[_idx] = val;
    }

  private:
    se::Node<FieldType> * _node; 
    int _idx; 
};


#endif
