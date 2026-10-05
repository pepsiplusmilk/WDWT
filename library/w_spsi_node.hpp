//
// Created by Alex Vladimirov on 05.10.2026.
//

#ifndef W_SPSI_NODE_H
#define W_SPSI_NODE_H
#include "weighted_spsi.hpp"
#include <array>
#include <cassert>
#include <vector>

template <class leaf_type, // implementation of leaf container
uint32_t B_leaf, // Number of element in leaf of B-tree is between B_leaf and 2B_leaf
uint_32_t B_fan_out> // Number of childrens of vertex is between B_fan_out + 1 and 2B_fan_out + 2
class w_spsi<leaf_type, B_leaf, B_fan_out>::node {
public:
  void free_mem() {
    if (has_leaves()) {
      for (uint32_t i = 0; i < nr_children_; ++i) delete leaves[i];

    } else {
      for (uint32_t i = 0; i < nr_children_; ++i) children[i]->free_mem();
      for (uint32_t i = 0; i < nr_children_; ++i) delete children[i];
    }
  }

private:
  static uint64_t free_capacity(const leaf_type& l) {
    assert(l.size() <= 2 * B_leaf);
    return 2 * B_leaf - l.size();
  }

  bool is_root() const {
    return parent_ == nullptr;
  }

  bool can_loose_child() const {
    return (is_root() || nr_children_ >= B_fan_out + 2);
  }

  bool has_leaves() const {
    return has_leaves_;
  }

  /*
   * helper functions for child search
   */
  uint64_t find_child(uint64_t i) const {
    uint64_t j = 0;
    while (subtree_sizes[j] <= i) {
      j++;
      assert(j < subtree_sizes.size());
    }
    return j;
  }

  uint64_t find_1(uint64_t x) const {
    uint64_t j = 0;
    while (!subtree_psums[j] || subtree_psums[j] < x) {
      j++;
      assert(j < subtree_psums.size());
    }
    return j;
  }

  uint64_t find_0(uint64_t x) const {
    uint64_t j = 0;
    while (subtree_sizes[j] - subtree_psums[j] < x) {
      j++;
      assert(j < subtree_psums.size());
    }
    return j;
  }

  void set_rank(uint32_t r) { rank_ = r; }
  void set_parent(node* new_p) {parent_ = new_p; }

  std::array<uint64_t, 2 * B_fan_out> subtree_sizes;

  std::vector<node*> children;
  std::vector<leaf_type*> leaves; // empty if vertex isn't leaf

  node* parent_ = nullptr;

  uint32_t rank_ = 0;
  uint32_t nr_children_ = 0;
  bool has_leaves_ = false;
};


#endif  // W_SPSI_NODE_H
