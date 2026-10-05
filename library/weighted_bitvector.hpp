//
// Created by Alex Vladimirov on 04.10.2026.
//

#ifndef WEIGHTED_BITVECTOR_H
#define WEIGHTED_BITVECTOR_H

#include "monoid.hpp"
#include <cstddef>

template <typename w_spsi, class monoid>
class weighted_bitvector {
  w_spsi spsi_tree;
  using w_t = typename monoid::value_type;
public:
  // wt primitives
  bool access(size_t i);
  size_t rank(bool b, size_t i) const;
  size_t select(bool b, size_t i) const;

  // aggregations
  w_t range_aggregate(size_t l, size_t r) const;
  void set_weight(size_t i, w_t w);

  // dyn
  void insert(size_t i, bool b, w_t w);
  void delete(size_t i);
};

#endif  // WEIGHTED_BITVECTOR_H
