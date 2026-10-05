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
uint32_t B_fan_out> // Number of childrens of vertex is between B_fan_out + 1 and 2B_fan_out + 2
class w_spsi<leaf_type, B_leaf, B_fan_out>::node {
public:
  using w_t = typename leaf_type::w_t;
  using monoid = typename leaf_type::monoid_type;
  static constexpr size_t MAX_CHILDREN = 2 * B_fan_out + 2;

  node() : nr_children_(1), has_leaves_(true) {
    leaves.push_back(new leaf_type());

    subtree_sizes[0] = 0;
    subtree_psums[0] = 0;
    subtree_aggs[0] = monoid::id;
  }

  node(std::vector<node*>&& c, node* P = nullptr, uint32_t rank = 0)
    : parent_(P), rank_(rank), nr_children_(c.size()), children(std::move(c)) {
    uint32_t r = 0;

    for (auto* child : children) {
      child->set_parent(this);
      child->set_rank(r++);
    }

    rebuild(0);
  }

  node(std::vector<leaf_type*>&& l, node* P = nullptr, uint32_t rank = 0)
    : parent_(P), rank_(rank), nr_children_(l.size()), has_leaves_(true), leaves(std::move(l)) {
    rebuild(0);
  }

  w_t range_aggregate(uint64_t l, uint64_t r) const {
    assert(l <= r && r < size());
    uint32_t jl = find_child(l);
    uint32_t jr = find_child(r);

    uint64_t prev_l = (jl == 0 ? 0 : subtree_sizes[jl - 1]);

    if (jl == jr) {
      if (has_leaves_) {
        return leaves[jl]->range_aggregate(l - prev_l, r - prev_l);
      }
      return children[jl]->range_aggregate(l - prev_l, r - prev_l);
    }

    uint64_t child_l_sz = (has_leaves_ ? leaves[jl]->size() : children[jl]->size());
    w_t res = (has_leaves_ ? leaves[jl]->range_aggregate(l - prev_l, child_l_sz - 1)
                           : children[jl]->range_aggregate(l - prev_l, child_l_sz - 1));

    for (uint32_t k = jl + 1; k < jr; ++k) {
      w_t mid_weight = (has_leaves_ ? leaves[k]->total_weight() : children[k]->total_weight());
      res = Monoid::combine(res, mid_weight);
    }

    // 3. Префикс крайнего правого ребенка
    uint64_t prev_r = subtree_sizes[jr - 1];
    w_t right_res = (has_leaves_ ? leaves[jr]->range_aggregate(0, r - prev_r)
                                 : children[jr]->range_aggregate(0, r - prev_r));

    return Monoid::combine(res, right_res);
  }

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
    assert(l.size() <= MAX_CHILDREN);
    return MAX_CHILDREN - l.size();
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

  void rebuild(size_t start_idx) {
    uint64_t sz = (start_idx == 0 ? 0 : subtree_sizes[start_idx - 1]);
    uint64_t psums = (start_idx == 0 ? 0 : subtree_psums[start_idx - 1]);
    w_t agg = (start_idx == 0 ? monoid::id : subtree_aggs[start_idx - 1]);

    for (size_t k = start_idx; k < nr_children_; ++k) {
      if (has_leaves_) {
        sz += leaves[k]->size();
        psums += leaves[k]->psum();

        agg = monoid::combine(agg, leaves[k]->total_weight());
      } else {
        sz += children[k]->size();
        psums += children[k]->total_ones();

        agg = monoid::combine(agg, children[k]->total_weight());
      }

      subtree_sizes[k] = sz;
      subtree_psums[k] = psums;
      subtree_aggs[k] = agg;
    }
  }

  void set_rank(uint32_t r) { rank_ = r; }
  void set_parent(node* new_p) {parent_ = new_p; }

  std::array<uint64_t, MAX_CHILDREN> subtree_sizes;

  std::array<uint64_t, MAX_CHILDREN> subtree_psums;
  std::array<w_t, MAX_CHILDREN> subtree_aggs;

  std::vector<node*> children;
  std::vector<leaf_type*> leaves; // empty if vertex isn't leaf

  node* parent_ = nullptr;

  uint32_t rank_ = 0;
  uint32_t nr_children_ = 0;
  bool has_leaves_ = false;
};


#endif  // W_SPSI_NODE_H
