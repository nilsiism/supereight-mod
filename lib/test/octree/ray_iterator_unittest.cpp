#include "utils/se_common.h"
#include "octree.hpp"
#include "gtest/gtest.h"
#include <vector>

typedef float testT;

template <>
struct voxel_traits<testT> {
  typedef float value_type;
  static inline value_type empty(){ return 0.f; }
  static inline value_type initValue(){ return 1.f; }
};
typedef Octree<testT> OctreeF;

class RayIteratorTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
      oct_.init(512, 5);

      p_   = Eigen::Vector3f(1.5f, 1.5f, 1.5f);
      dir_ = Eigen::Vector3f(0.5, 0.5, 0.5).normalized();

      // ensure stepsize is big enough to get distinct blocks
      const float stepsize = 2 * (oct_.dim()/oct_.size() * OctreeF::blockSide);
      const float voxelsize = oct_.dim()/oct_.size();

      const int num_blocks = 4;
      float t = 0.6f;
      for(int i = 0; i < num_blocks; ++i, t += stepsize) {
        const Eigen::Vector3f tmp = p_ + t * dir_;
        const Eigen::Vector3i vox = ((p_ + t * dir_)/voxelsize).cast<int> ();

        // hash to VoxelBlocks
        se::key_t key = oct_.hash(vox(0), vox(1), vox(2)) ;
        alloc_list_.push_back(key);
      }

      oct_.allocate(alloc_list_.data(), alloc_list_.size());
    }
  OctreeF oct_;
  Eigen::Vector3f p_;
  Eigen::Vector3f dir_;
  std::vector<se::key_t> alloc_list_;
};

TEST_F(RayIteratorTest, FetchAlongRay) {
  se::ray_iterator<OctreeF::value_type> it(oct_, p_, dir_, 0.4, 4.0f); 
  int i = 0;
  se::VoxelBlock<testT> * current;
  while(current = it.next()) {
    ASSERT_LT(i, alloc_list_.size());
    ASSERT_EQ(current->code, alloc_list_[i]);
    i++; 
  }
  ASSERT_EQ(i, alloc_list_.size());
}
