#ifndef OCTREE_COLLISION_HPP
#define OCTREE_COLLISION_HPP
#include "node.hpp"
#include "octree.hpp"
#include "geometry/aabb_collision.hpp"

/*! \brief Perform a collision test betweenb the input octree map and the 
 * input axis aligned bounding box bbox of extension side. The test function 
 * test takes as input a voxel value and returns a collision_status. This is
 * used to distinguish between seen-empty voxels and occupied voxels.
 * \param map octree map
 * \param bbox test bounding box lower bottom corner 
 * \param side extension in number of voxels of the bounding box
 * \param test function that takes a voxel and returns a collision_status value
 */

template <typename FieldType, typename TestVoxelF>
collision_status collision_test(const Octree<FieldType>& map, 
    const int3 bbox, const int3 side, TestVoxelF test) {

  typedef struct stack_entry { 
    Node<FieldType>* node_ptr;
    int3 coordinates;
    int side;
  } stack_entry;

  constexpr int3 offsets[8] = {{0, 0 ,0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
    {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};

  stack_entry stack[Octree<FieldType>::max_depth*8 + 1];
  size_t stack_idx = 0;

  Node<FieldType>* n = map.root();
  if(!n) return collision_status::unseen;

  stack_entry current;
  current.node_ptr = n;
  current.side = map.size();
  current.coordinates = {0, 0, 0};
  stack[stack_idx++] = current;
  collision_status status = collision_status::empty;

  while(stack_idx != 0){
    Node<FieldType>* node = current.node_ptr;

    if(node->isLeaf()){
      return collision_status::occupied;
    } 

    status = update_status(status, test(node->value_));
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
        stack[stack_idx++] = child_descr;
      }
    }
    current = stack[--stack_idx]; 
  }
  return status;
}
#endif
