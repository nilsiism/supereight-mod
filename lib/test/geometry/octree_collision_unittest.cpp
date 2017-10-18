#include "math_utils.h"
#include "geometry/octree_collision.hpp"
#include "geometry/aabb_collision.hpp"
#include "algorithms/mapping.hpp"
#include "utils/morton_utils.hpp"
#include "octree.hpp"
#include "gtest/gtest.h"

typedef float testT;

template <>
struct voxel_traits<testT> {
  typedef float ComputeType;
  typedef float StoredType;
  static inline ComputeType empty(){ return 0.f; }
  static inline ComputeType initValue(){ return 1.f; }
  static inline StoredType translate(const ComputeType value) {
     return value;
  }
};

collision_status test_voxel(const voxel_traits<testT>::ComputeType & val) {
  if(val == voxel_traits<testT>::initValue()) return collision_status::unseen;
  if(val == 10.f) return collision_status::empty;
  return collision_status::occupied;
};

class OctreeCollisionTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      oct_.init(256, 5);
      const int3 blocks[1] = {{56, 12, 254}};
      unsigned int alloc_list[1];
      alloc_list[0] = oct_.hash(blocks[0].x, blocks[0].y, blocks[0].z);

      auto update = [](Node<testT> * n){
        uint3 coords = unpack_morton(n->code);
        /* Empty for coords above the below values, 
         * except where leaves are allocated.
         */
        if(coords.x >= 48 && coords.y >= 0 && coords.z >= 240) {
          n->value_ = 10.f;
        }
      };
      oct_.alloc_update(alloc_list, 1, 6, update);
      algorithms::integratePass(oct_.getNodesBuffer(), oct_.getNodesBuffer().size(),
          update);
    }

  typedef Octree<testT> OctreeF;
  OctreeF oct_;
};

TEST_F(OctreeCollisionTest, TotallyUnseen) {

  leaf_iterator<testT> it(oct_);
  typedef std::tuple<int3, int, typename Octree<testT>::compute_type> it_result;
  it_result node = it.next();
  for(int i = 128; std::get<1>(node) > 0; node = it.next(), i /= 2){
    const int3 coords = std::get<0>(node);
    const int side = std::get<1>(node);
    const Octree<testT>::compute_type val = std::get<2>(node);
    printf("Node's coordinates: (%d, %d, %d), side %d, value %.2f\n", 
        coords.x, coords.y, coords.z, side, val);
    EXPECT_EQ(side, i);
  }

  const int3 test_bbox = {23, 0, 100};
  const int3 width = {2, 2, 2};

  const collision_status collides = collision_test(oct_, test_bbox, width, 
      test_voxel);
  ASSERT_EQ(collides, collision_status::unseen);
}

TEST_F(OctreeCollisionTest, PartiallyUnseen) {
  const int3 test_bbox = {47, 0, 239};
  const int3 width = {6, 6, 6};
  const collision_status collides = collision_test(oct_, test_bbox, width, 
      test_voxel);
  ASSERT_EQ(collides, collision_status::unseen);
}

TEST_F(OctreeCollisionTest, Empty) {
  const int3 test_bbox = {49, 1, 242};
  const int3 width = {1, 1, 1};
  const collision_status collides = collision_test(oct_, test_bbox, width, 
      test_voxel);
  ASSERT_EQ(collides, collision_status::empty);
}

// TEST_F(OctreeCollisionTest, CollisionPossible){
//   const int3 test_bbox = {54, 10, 249};
//   const int3 width = {5, 5, 3};
// 
//   const collision_status collides = collision_test(oct_, test_bbox, width, 
//       test_voxel);
//   ASSERT_EQ(collides, collision_status::occupied);
// }
