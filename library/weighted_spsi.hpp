//
// Modified by Alex Vladimirov on 04.10.2026.
//

// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * spsi.hpp
 *
 *  Created on: Oct 19, 2015
 *      Author: nico
 *
 *  Searchable partial sums with insert.
 *
 *  represents a vector of integers I_0, ..., I_(n-1) >= 0
 *  supports random access, set, partial sum, search, insert, update
 * (increment/decrement).
 *
 *  The structure is a B+-tree. This improves data locality and space
 * efficiency.
 *
 */

#ifndef W_SPSI_H
#define W_SPSI_H

#include <cstdint>
#include <cassert>

template <class leaf_type, // implementation of leaf container
uint32_t B_leaf, // Number of element in leaf of B-tree is between B_leaf and 2B_leaf
uint32_t B_fan_out> // Number of childrens of vertex is between B_fan_out + 1 and 2B_fan_out + 2
class w_spsi {
public:
  using w_t = typename leaf_type::w_t;
  using monoid = typename leaf_type::monoid_type;

  w_spsi();
  ~w_spsi();

  w_spsi(const w_spsi&) = delete;
  w_spsi& operator=(const w_spsi&) = delete;
  w_spsi(w_spsi&& other) noexcept : root_(other.root_) { other.root_ = nullptr; }
  w_spsi& operator=(w_spsi&& other) noexcept;

  uint64_t size() const { return root_ ? root_->size() : 0; }
  bool empty() const { return size() == 0; }
  uint64_t total_ones() const { return root_ ? root_->total_ones() : 0; }
  w_t total_weight() const { return root_ ? root_->total_weight() : monoid::id; }

  bool at(uint64_t i) const { assert(root_); return root_->at(i); }
  w_t weight_at(uint64_t i) const { assert(root_); return root_->weight_at(i); }

  uint64_t rank_1(uint64_t i) const { assert(root_); return root_->rank_1(i); }
  uint64_t rank_0(uint64_t i) const { assert(root_); return root_->rank_0(i); }
  uint64_t select_1(uint64_t x) const { assert(root_); return root_->select_1(x); }
  uint64_t select_0(uint64_t x) const { assert(root_); return root_->select_0(x); }

  w_t range_aggregate(uint64_t l, uint64_t r) const {
    if (l > r || r >= size() || !root_) return monoid::id;
    return root_->range_aggregate(l, r);
  }

  void set_weight(uint64_t i, w_t w) {
    assert(root_);
    root_->set_weight(i, w);
  }

  void range_set_weight(uint64_t l, uint64_t r, w_t w) {
    if (l > r || r >= size() || !root_) return;
    root_->range_set_weight(l, r, w);
  }

  void insert(uint64_t i, bool bit, w_t weight) {
    assert(root_);

    node* new_root = root_->insert(i, bit, weight);

    if (new_root) {
      root_ = new_root;
    }
  }

  void push_back(bool bit, w_t weight) {
    insert(size(), bit, weight);
  }

  void remove(uint64_t i) {
    assert(root_);

    node* new_root = root_->remove(i);

    if (new_root) {
      delete root_;
      root_ = new_root;
    }
  }

  size_t size_in_bits() const {
    size_t bits = sizeof(*this) * 8;
    if (root_) {
      bits += root_->size_in_bits();
    }
    return bits;
  }

  size_t size_in_bytes() const {
    return (size_in_bits() + 7) / 8;
  }

private:
  class node;
  node* root_ = nullptr; // B-tree root
};

#include "w_spsi_node.hpp"

template <class leaf_type, uint32_t B_leaf, uint32_t B_fan_out>
w_spsi<leaf_type, B_leaf, B_fan_out>::w_spsi() : root_(new node()) {
}

template <class leaf_type, uint32_t B_leaf, uint32_t B_fan_out>
w_spsi<leaf_type, B_leaf, B_fan_out>::~w_spsi() {
  if (root_) {
    root_->free_mem();
    delete root_;
    root_ = nullptr;
  }
}

template <class leaf_type, uint32_t B_leaf, uint32_t B_fan_out>
w_spsi<leaf_type, B_leaf, B_fan_out>& w_spsi<leaf_type, B_leaf, B_fan_out>::operator=(w_spsi&& other) noexcept {
  if (this != &other) {
    if (root_) {
      root_->free_mem();
      delete root_;
    }
    root_ = other.root_;
    other.root_ = nullptr;
  }
  return *this;
}

#endif // W_SPSI_H