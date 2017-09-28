#include "octree.hpp"
#include "math_utils.h"
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

class MultiscaleTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      oct_.init(512, 5);

    }

  typedef Octree<testT> OctreeF;
  OctreeF oct_;
};

TEST_F(MultiscaleTest, Init) {
  EXPECT_EQ(oct_.get(137, 138, 130), voxel_traits<testT>::initValue());
}

TEST_F(MultiscaleTest, PlainAlloc) {
  const int3 blocks[2] = {{56, 12, 254}, {87, 32, 423}};
  unsigned int alloc_list[2];
  for(int i = 0; i < 2; ++i) {
    alloc_list[i] = oct_.hash(blocks[i].x, blocks[i].y, blocks[i].z);
  }
  oct_.allocate(alloc_list, 2);

  oct_.set(56, 12, 254, 3.f);

  EXPECT_EQ(oct_.get(56, 12, 254), 3.f);
  EXPECT_EQ(oct_.get(106, 12, 254), voxel_traits<testT>::initValue());
  EXPECT_NE(oct_.get(106, 12, 254), 3.f);
}

TEST_F(MultiscaleTest, ScaledAlloc) {
  const int3 blocks[2] = {{56, 12, 254}, {87, 32, 423}};
  unsigned int alloc_list[2];
  for(int i = 0; i < 2; ++i) {
    alloc_list[i] = oct_.hash(blocks[i].x, blocks[i].y, blocks[i].z);
  }

  auto update = [](Node<testT> * n, const int , const int , const int ){
    n->value_ = 10.f;
  };
  oct_.alloc_update(alloc_list, 2, 5, update);
  Node<testT>* n = oct_.fetch_octant(56, 12, 254, 5);
  ASSERT_TRUE(n != NULL);
  n->value_ = 10.f;
  EXPECT_EQ(oct_.get(56, 12, 254), 10.f);
}
