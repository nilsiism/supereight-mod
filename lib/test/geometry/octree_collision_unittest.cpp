#include "math_utils.h"
#include "geometry/aabb_collision.hpp"
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

class OctreeCollisionTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      oct_.init(512, 5);
      const int3 blocks[10] = {{56, 12, 254}, {87, 32, 423}, {128, 128, 128},
      {136, 128, 128}, {128, 136, 128}, {136, 136, 128}, 
      {128, 128, 136}, {136, 128, 136}, {128, 136, 136}, {136, 136, 136}};
      unsigned int alloc_list[10];
      for(int i = 0; i < 10; ++i) {
        alloc_list[i] = oct_.hash(blocks[i].x, blocks[i].y, blocks[i].z);
      }
      oct_.allocate(alloc_list, 10);
    }

  typedef Octree<testT> OctreeF;
  OctreeF oct_;
};

TEST_F(OctreeCollisionTest, CollisionMiss){
  const int3 test_bbox = {100, 100, 100};
  const int3 width = {5, 5, 5};

  const int collides = oct_.collision_test(test_bbox, width);
  ASSERT_EQ(collides, 0);
}

TEST_F(OctreeCollisionTest, CollisionPossible){
  const int3 test_bbox = {54, 10, 256};
  const int3 width = {5, 5, 5};

  const int collides = oct_.collision_test(test_bbox, width);
  ASSERT_EQ(collides, 1);
}
