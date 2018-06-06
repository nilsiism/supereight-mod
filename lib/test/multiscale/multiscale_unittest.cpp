#include "octree.hpp"
#include "utils/se_common.h"
#include "gtest/gtest.h"

typedef float testT;

template <>
struct voxel_traits<testT> {
  typedef float value_type;
  static inline value_type empty(){ return 0.f; }
  static inline value_type initValue(){ return 1.f; }
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
  const Eigen::Vector3i blocks[2] = {{56, 12, 254}, {87, 32, 423}};
  se::key_t alloc_list[2];
  for(int i = 0; i < 2; ++i) {
    alloc_list[i] = oct_.hash(blocks[i](0), blocks[i](1), blocks[i](2));
  }
  oct_.allocate(alloc_list, 2);

  oct_.set(56, 12, 254, 3.f);

  EXPECT_EQ(oct_.get(56, 12, 254), 3.f);
  EXPECT_EQ(oct_.get(106, 12, 254), voxel_traits<testT>::initValue());
  EXPECT_NE(oct_.get(106, 12, 254), 3.f);
}

TEST_F(MultiscaleTest, ScaledAlloc) {
  const Eigen::Vector3i blocks[2] = {{200, 12, 25}, {87, 32, 423}};
  se::key_t alloc_list[2];
  for(int i = 0; i < 2; ++i) {
    alloc_list[i] = oct_.hash(blocks[i](0), blocks[i](1), blocks[i](2), 5);
  }

  oct_.allocate(alloc_list, 2);
  Node<testT>* n = oct_.fetch_octant(87, 32, 420, 5);
  ASSERT_TRUE(n != NULL);
  n->value_[0] = 10.f;
  EXPECT_EQ(oct_.get(87, 32, 420), 10.f);
}

TEST_F(MultiscaleTest, Iterator) {
  const Eigen::Vector3i blocks[1] = {{56, 12, 254}};
  se::key_t alloc_list[1];
  alloc_list[0] = oct_.hash(blocks[0](0), blocks[0](1), blocks[0](2));

  oct_.allocate(alloc_list, 1);
  leaf_iterator<testT> it(oct_);

  typedef std::tuple<Eigen::Vector3i, int, typename Octree<testT>::value_type> it_result;
  it_result node = it.next();
  for(int i = 256; std::get<1>(node) > 0; node = it.next(), i /= 2){
    const Eigen::Vector3i coords = std::get<0>(node);
    const int side = std::get<1>(node);
    const Octree<testT>::value_type val = std::get<2>(node);
    EXPECT_EQ(side, i);
  }
}

TEST_F(MultiscaleTest, ChildrenMaskTest) {
  const Eigen::Vector3i blocks[10] = {{56, 12, 254}, {87, 32, 423}, {128, 128, 128},
    {136, 128, 128}, {128, 136, 128}, {136, 136, 128}, 
    {128, 128, 136}, {136, 128, 136}, {128, 136, 136}, {136, 136, 136}};
  se::key_t alloc_list[10];
  for(int i = 0; i < 10; ++i) {
    alloc_list[i] = oct_.hash(blocks[i](0), blocks[i](1), blocks[i](2), 5);
  }

  oct_.allocate(alloc_list, 10);
  const MemoryPool<Node<testT> >& nodes = oct_.getNodesBuffer();
  const size_t num_nodes = nodes.size();
  for(size_t i = 0; i < num_nodes; ++i) {
    Node<testT>* n = nodes[i];
    for(int c = 0; c < 8; ++c) {
      if(n->child(c)) {
        ASSERT_TRUE(n->children_mask_ & (1 << c));
      }
    }
  } 
}

TEST_F(MultiscaleTest, OctantAlloc) {
  const Eigen::Vector3i blocks[10] = {{56, 12, 254}, {87, 32, 423}, {128, 128, 128},
    {136, 128, 128}, {128, 136, 128}, {136, 136, 128}, 
    {128, 128, 136}, {136, 128, 136}, {128, 136, 136}, {136, 136, 136}};
  se::key_t alloc_list[10];
  for(int i = 0; i < 10; ++i) {
    alloc_list[i] = oct_.hash(blocks[i](0), blocks[i](1), blocks[i](2));
  }

  alloc_list[2] = alloc_list[2] | 3;
  alloc_list[9] = alloc_list[2] | 5;
  oct_.allocate(alloc_list, 10);
  Node<testT> * octant = oct_.fetch_octant(blocks[4](0), blocks[4](1),
      blocks[4](2), 3);
  ASSERT_TRUE(octant != NULL);
  octant = oct_.fetch_octant(blocks[9](0), blocks[9](1),
      blocks[9](2), 6);
  ASSERT_TRUE(octant == NULL);
}
