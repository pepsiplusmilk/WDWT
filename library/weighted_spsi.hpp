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

template <class leaf_type, // implementation of leaf container
uint32_t B_leaf, // Number of element in leaf of B-tree is between B_leaf and 2B_leaf
uint_32_t B_fan_out> // Number of childrens of vertex is between B_fan_out + 1 and 2B_fan_out + 2
class w_spsi {
  class node;
  node* root = nullptr; // B-tree root

public:
  w_spsi();
};

#endif //W_SPSI_H
