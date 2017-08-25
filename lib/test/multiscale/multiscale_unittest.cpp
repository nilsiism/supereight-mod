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

class InterpTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      oct_.init(512, 5);
    }

  typedef Octree<testT> OctreeF;
  OctreeF oct_;
};

TEST_F(InterpTest, Init) {
  EXPECT_EQ(oct_.get(137, 138, 130), voxel_traits<testT>::empty());
}
