//
// Created by Alex Vladimirov on 04.10.2026.
//

#ifndef WEIGHTED_BITVECTOR_H
#define WEIGHTED_BITVECTOR_H

#include "monoid.hpp"
#include <cstddef>
#include <cstdint>

template <typename w_spsi, class monoid>
class weighted_bitvector {
  w_spsi spsi_tree;
  using w_t = typename monoid::value_type;
public:
  weighted_bitvector() = default;

  size_t size() const { return spsi_tree.size(); }
  bool empty() const { return spsi_tree.empty(); }

  // wt primitives
  bool access(size_t i) const { return spsi_tree.at(i); }
  w_t weight_at(size_t i) const { return spsi_tree.weight_at(i); }

  size_t rank(bool b, size_t i) const {
    if (spsi_tree.size() == 0 || i == static_cast<size_t>(-1)) return 0;
    if (i >= spsi_tree.size()) i = spsi_tree.size() - 1;
    return b ? spsi_tree.rank_1(i) : spsi_tree.rank_0(i);
  }

  size_t select(bool b, size_t i) const {
    return b ? spsi_tree.select_1(i) : spsi_tree.select_0(i);
  }

  // aggregations
  w_t range_aggregate(size_t l, size_t r) const {
    if (l > r || r >= spsi_tree.size()) return monoid::id;
    return spsi_tree.range_aggregate(l, r);
  }

  void set_weight(size_t i, w_t w) {
    spsi_tree.set_weight(i, w);
  }

  void range_set_weight(size_t l, size_t r, w_t w) {
    if (l <= r && r < spsi_tree.size()) {
      spsi_tree.range_set_weight(l, r, w);
    }
  }

  // dyn
  void insert(size_t i, bool b, w_t w) {
    spsi_tree.insert(i, b, w);
  }

  void push_back(bool b, w_t w) {
    spsi_tree.push_back(b, w);
  }

  void remove(size_t i) {
    spsi_tree.remove(i);
  }
};

#endif // WEIGHTED_BITVECTOR_H