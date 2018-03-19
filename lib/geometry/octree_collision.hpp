#ifndef OCTREE_COLLISION_HPP
#define OCTREE_COLLISION_HPP
#include "node.hpp"
#include "octree.hpp"
#include "geometry/aabb_collision.hpp"

enum class collision_status {
  occupied,
  unseen,
  empty
};

/*! \brief Implements a simple state machine to update the collision status.
 * The importance order is given as follows in ascending order: 
 * Empty, Unseen, Occupied.
 * \param previous_status
 * \param new_status 
 */
inline collision_status update_status(const collision_status previous_status, 
    const collision_status new_status) {
  switch(previous_status) {
    case collision_status::unseen:
      if(new_status != collision_status::occupied)
        return previous_status;
      else 
        return new_status;
      break;
    case collision_status::occupied:
      return previous_status;
      break;
    default:
      return new_status;
      break;
  }
}

/*! \brief Perform a collision test for each voxel value in the input voxel 
 * block. The test function test takes as input a voxel value and returns a
 * collision_status. This is used to distinguish between seen-empty voxels and
 * occupied voxels.
 * \param block voxel block of type FieldType
 * \param test function that takes a voxel and returns a collision_status value
 */
template <typename FieldType, typename TestVoxelF>
collision_status collides_with(const VoxelBlock<FieldType>* block, 
    const int3 bbox, const int3 side, TestVoxelF test) {
  collision_status status = collision_status::empty;
  const int3 blockCoord = block->coordinates();
  int x, y, z, blockSide; 
  blockSide = (int) VoxelBlock<FieldType>::side;
  int xlast = blockCoord.x + blockSide;
  int ylast = blockCoord.y + blockSide;
  int zlast = blockCoord.z + blockSide;
  for(z = blockCoord.z; z < zlast; ++z){
    for (y = blockCoord.y; y < ylast; ++y){
      for (x = blockCoord.x; x < xlast; ++x){

        typename VoxelBlock<FieldType>::compute_type value;
        const int3 vox = make_int3(x, y, z);
        if(!geometry::aabb_aabb_collision(bbox, side, 
          vox, make_int3(1))) continue;
        value = block->data(make_int3(x, y, z));
        status = update_status(status, test(value));
      }
    }
  }
  return status;
}

/*! \brief Perform a collision test between the input octree map and the 
 * input axis aligned bounding box bbox of extension side. The test function 
 * test takes as input a voxel value and returns a collision_status. This is
 * used to distinguish between seen-empty voxels and occupied voxels.
 * \param map octree map
 * \param bbox test bounding box lower bottom corner 
 * \param side extension in number of voxels of the bounding box
 * \param test function that takes a voxel and returns a collision_status value
 */

template <typename FieldType, typename TestVoxelF>
collision_status collides_with(const Octree<FieldType>& map, 
    const int3 bbox, const int3 side, TestVoxelF test) {

  typedef struct stack_entry { 
    Node<FieldType>* node_ptr;
    int3 coordinates;
    int side;
    typename Node<FieldType>::compute_type parent_val;
  } stack_entry;

  constexpr int3 offsets[8] = {{0, 0 ,0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
    {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};

  stack_entry stack[Octree<FieldType>::max_depth*8 + 1];
  size_t stack_idx = 0;

  Node<FieldType>* node = map.root();
  if(!node) return collision_status::unseen;

  stack_entry current;
  current.node_ptr = node;
  current.side = map.size();
  current.coordinates = {0, 0, 0};
  stack[stack_idx++] = current;
  collision_status status = collision_status::empty;

  while(stack_idx != 0){
    node = current.node_ptr;

    if(node->isLeaf()){
      status = collides_with(static_cast<VoxelBlock<FieldType>*>(node), 
          bbox, side, test);
    } 

    if(node->children_mask_ == 0) {
       current = stack[--stack_idx]; 
       continue;
    }

    for(int i = 0; i < 8; ++i){
      Node<FieldType>* child = node->child(i);
      stack_entry child_descr;
      child_descr.node_ptr = NULL;
      child_descr.side = current.side / 2;
      child_descr.coordinates = 
          make_int3(current.coordinates.x + child_descr.side*offsets[i].x,
                    current.coordinates.y + child_descr.side*offsets[i].y,
                    current.coordinates.z + child_descr.side*offsets[i].z);
      const bool overlaps = geometry::aabb_aabb_collision(bbox, side, 
          child_descr.coordinates, 
              make_int3(child_descr.side));

      if(overlaps && child != NULL) {
        child_descr.node_ptr = child;
        child_descr.parent_val = node->value_[0];
        stack[stack_idx++] = child_descr;
      } else if(overlaps && child == NULL) {
        status = update_status(status, test(node->value_[0]));
      }
    }
    current = stack[--stack_idx]; 
  }
  return status;
}
#endif
