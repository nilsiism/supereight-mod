#include "utils/se_common.h"
#include "octree.hpp"
#include "gtest/gtest.h"
#include <vector>

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
        octlib::key_t key = oct_.hash(vox(0), vox(1), vox(2)) ;
        alloc_list_.push_back(key);
      }

      oct_.alloc_update(alloc_list_.data(), alloc_list_.size());
    }
  OctreeF oct_;
  Eigen::Vector3f p_;
  Eigen::Vector3f dir_;
  std::vector<octlib::key_t> alloc_list_;
};

TEST_F(RayIteratorTest, FetchAlongRay) {
  octlib::ray_iterator<OctreeF::compute_type> it(oct_, p_, dir_, 0.4, 4.0f); 
  int i = 0;
  VoxelBlock<testT> * current;
  while(current = it.next()) {
    ASSERT_LT(i, alloc_list_.size());
    ASSERT_EQ(current->code, alloc_list_[i]);
    i++; 
  }
  ASSERT_EQ(i, alloc_list_.size());
}

TEST_F(RayIteratorTest, CheckTminTmax) {
  // const float blocksize = OctreeF::blockSide * (oct_.dim()/oct_.size());
  // const float diag = sqrtf((blocksize*blocksize) + (blocksize*blocksize));
  // octlib::ray_iterator<OctreeF::compute_type> it(oct_, p_, dir_, 0.4, 4.0f); 

  // int i = 0;
  // VoxelBlock<testT> * current;
  // while(current = it.next()) {
  //   float diff = it.tcmax() - it.tcmin();
  //   std::cout << diff << std::endl;
  //   ASSERT_LT(diff, diag);
  //   i++; 
  // }
  // ASSERT_EQ(i, alloc_list_.size());
}
