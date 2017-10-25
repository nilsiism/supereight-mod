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
    alloc_list[i] = oct_.hash(blocks[i].x, blocks[i].y, blocks[i].z, 5);
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

TEST_F(MultiscaleTest, Iterator) {
  const int3 blocks[1] = {{56, 12, 254}};
  unsigned int alloc_list[1];
  alloc_list[0] = oct_.hash(blocks[0].x, blocks[0].y, blocks[0].z);

  auto update = [](Node<testT> * n, const int , const int , const int ){
    n->value_ = 10.f;
  };
  oct_.alloc_update(alloc_list, 1, 7, update);
  leaf_iterator<testT> it(oct_);

  typedef std::tuple<int3, int, typename Octree<testT>::compute_type> it_result;
  it_result node = it.next();
  for(int i = 256; std::get<1>(node) > 0; node = it.next(), i /= 2){
    const int3 coords = std::get<0>(node);
    const int side = std::get<1>(node);
    const Octree<testT>::compute_type val = std::get<2>(node);
    EXPECT_EQ(side, i);
  }
}

TEST_F(MultiscaleTest, ChildrenMaskTest) {
  const int3 blocks[10] = {{56, 12, 254}, {87, 32, 423}, {128, 128, 128},
    {136, 128, 128}, {128, 136, 128}, {136, 136, 128}, 
    {128, 128, 136}, {136, 128, 136}, {128, 136, 136}, {136, 136, 136}};
  unsigned int alloc_list[10];
  for(int i = 0; i < 10; ++i) {
    alloc_list[i] = oct_.hash(blocks[i].x, blocks[i].y, blocks[i].z, 5);
  }

  auto update = [](Node<testT> * n, const int , const int , const int ){
    n->value_ = 10.f;
  };
  oct_.alloc_update(alloc_list, 10, 5, update);
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
  const int3 blocks[10] = {{56, 12, 254}, {87, 32, 423}, {128, 128, 128},
    {136, 128, 128}, {128, 136, 128}, {136, 136, 128}, 
    {128, 128, 136}, {136, 128, 136}, {128, 136, 136}, {136, 136, 136}};
  unsigned int alloc_list[10];
  for(int i = 0; i < 10; ++i) {
    alloc_list[i] = oct_.hash(blocks[i].x, blocks[i].y, blocks[i].z);
  }

  alloc_list[2] = alloc_list[2] | 3;
  alloc_list[9] = alloc_list[2] | 5;
  auto update = [](Node<testT> * n, const int , const int , const int ){
    n->value_ = 10.f;
  };
  oct_.alloc_update(alloc_list, 10, 7, update);
  Node<testT> * octant = oct_.fetch_octant(blocks[4].x, blocks[4].y,
      blocks[4].z, 3);
  ASSERT_TRUE(octant != NULL);
  octant = oct_.fetch_octant(blocks[9].x, blocks[9].y,
      blocks[9].z, 6);
  ASSERT_TRUE(octant == NULL);
}
