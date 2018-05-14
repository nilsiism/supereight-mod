#include "utils/eigen_helper.h"
#include "geometry/octree_collision.hpp"
#include "geometry/aabb_collision.hpp"
#include "utils/morton_utils.hpp"
#include "octree.hpp"
#include "functors/axis_aligned_functor.hpp"
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
      const Eigen::Vector3i blocks[1] = {{56, 12, 254}};
      octlib::key_t alloc_list[1];
      alloc_list[0] = oct_.hash(blocks[0](0), blocks[0](1), blocks[0](2));
      oct_.alloc_update(alloc_list, 1);

      auto set_to_ten = [](auto& handler, const Eigen::Vector3i& coords) {
        if((coords.array() >= Eigen::Vector3i(48, 0, 240).array()).all()){
          handler.set(10.f);
        }
      };
      iterators::functor::axis_aligned<testT, Octree, decltype(set_to_ten)> 
        funct(oct_, set_to_ten);
      funct.apply();
    }

  typedef Octree<testT> OctreeF;
  OctreeF oct_;
};

TEST_F(OctreeCollisionTest, TotallyUnseen) {

  leaf_iterator<testT> it(oct_);
  typedef std::tuple<Eigen::Vector3i, int, typename Octree<testT>::compute_type> it_result;
  it_result node = it.next();
  for(int i = 128; std::get<1>(node) > 0; node = it.next(), i /= 2){
    const Eigen::Vector3i coords = std::get<0>(node);
    const int side = std::get<1>(node);
    const Octree<testT>::compute_type val = std::get<2>(node);
    printf("Node's coordinates: (%d, %d, %d), side %d, value %.2f\n", 
        coords(0), coords(1), coords(2), side, val);
    EXPECT_EQ(side, i);
  }

  const Eigen::Vector3i test_bbox = {23, 0, 100};
  const Eigen::Vector3i width = {2, 2, 2};

  const collision_status collides = collides_with(oct_, test_bbox, width, 
      test_voxel);
  ASSERT_EQ(collides, collision_status::unseen);
}

TEST_F(OctreeCollisionTest, PartiallyUnseen) {
  const Eigen::Vector3i test_bbox = {47, 0, 239};
  const Eigen::Vector3i width = {6, 6, 6};
  const collision_status collides = collides_with(oct_, test_bbox, width, 
      test_voxel);
  ASSERT_EQ(collides, collision_status::unseen);
}

TEST_F(OctreeCollisionTest, Empty) {
  const Eigen::Vector3i test_bbox = {49, 1, 242};
  const Eigen::Vector3i width = {1, 1, 1};
  const collision_status collides = collides_with(oct_, test_bbox, width, 
      test_voxel);
  ASSERT_EQ(collides, collision_status::empty);
}

TEST_F(OctreeCollisionTest, Collision){
  const Eigen::Vector3i test_bbox = {54, 10, 249};
  const Eigen::Vector3i width = {5, 5, 3};

  auto update = [](auto& handler, const Eigen::Vector3i& coords) {
      handler.set(2.f);
  };
  iterators::functor::axis_aligned<testT, Octree, decltype(update)> 
    funct(oct_, update);
  funct.apply();
 
  const collision_status collides = collides_with(oct_, test_bbox, width, 
      test_voxel);
  ASSERT_EQ(collides, collision_status::occupied);
}

TEST_F(OctreeCollisionTest, CollisionFreeLeaf){
  // Allocated block: {56, 8, 248};
  const Eigen::Vector3i test_bbox = {61, 13, 253};
  const Eigen::Vector3i width = {2, 2, 2};

  /* Update leaves as occupied node */
  VoxelBlock<testT> * block = oct_.fetch(56, 12, 254);
  const Eigen::Vector3i blockCoord = block->coordinates();
  int x, y, z, blockSide; 
  blockSide = (int) VoxelBlock<testT>::side;
  int xlast = blockCoord(0) + blockSide;
  int ylast = blockCoord(1) + blockSide;
  int zlast = blockCoord(2) + blockSide;
  for(z = blockCoord(2); z < zlast; ++z){
    for (y = blockCoord(1); y < ylast; ++y){
      for (x = blockCoord(0); x < xlast; ++x){
        if(x < xlast/2 && y < ylast/2 && z < zlast/2)
          block->data(Eigen::Vector3i(x, y, z), 2.f);
        else
          block->data(Eigen::Vector3i(x, y, z), 10.f);

      }
    }
  }

  const collision_status collides = collides_with(oct_, test_bbox, width, 
      test_voxel);
  ASSERT_EQ(collides, collision_status::empty);
}
