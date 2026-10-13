//
// Created by Alex Vladimirov on 04.10.2026.
//

#ifndef WEIGHTED_WT_STRING_H
#define WEIGHTED_WT_STRING_H

#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <vector>
#include "monoid.hpp"
#include "w_pbv.hpp"
#include "weighted_spsi.hpp"
#include "weighted_bitvector.hpp"

template <typename w_t, size_t CHUNK_SIZE = 1024>
class base_leaf_storage {
public:
  struct Chunk {
    std::vector<char> chars;
    std::vector<w_t> weights;

    Chunk() {
      chars.reserve(CHUNK_SIZE * 2);
      weights.reserve(CHUNK_SIZE * 2);
    }

    size_t size() const {
      return chars.size();
    }

    bool empty() const {
      return size() == 0;
    }

    size_t size_in_bits() const {
      size_t bits = sizeof(*this) * 8;
      bits += chars.capacity() * sizeof(char) * 8;
      bits += weights.capacity() * sizeof(w_t) * 8;
      return bits;
    }
  };

  std::vector<Chunk> chunks;
  size_t total_size = 0;

  base_leaf_storage() {
    chunks.emplace_back();
  }

  size_t size() const {
    return total_size;
  }

  void push_back(char c, w_t w) {
    if (chunks.empty()) chunks.emplace_back();
    if (chunks.back().size() >= 2 * CHUNK_SIZE) {
      chunks.emplace_back();
    }
    chunks.back().chars.push_back(c);
    chunks.back().weights.push_back(w);
    total_size++;
  }

  void insert(size_t i, char c, w_t w) {
    assert(i <= total_size);
    if (i == total_size) {
      push_back(c, w);
      return;
    }
    size_t cur = 0;
    for (size_t b_idx = 0; b_idx < chunks.size(); ++b_idx) {
      if (i <= cur + chunks[b_idx].size()) {
        size_t offset = i - cur;
        chunks[b_idx].chars.insert(chunks[b_idx].chars.begin() + offset, c);
        chunks[b_idx].weights.insert(chunks[b_idx].weights.begin() + offset, w);
        total_size++;

        if (chunks[b_idx].size() >= 2 * CHUNK_SIZE) {
          Chunk new_chunk;
          size_t mid = chunks[b_idx].size() / 2;
          new_chunk.chars.assign(chunks[b_idx].chars.begin() + mid, chunks[b_idx].chars.end());
          new_chunk.weights.assign(chunks[b_idx].weights.begin() + mid, chunks[b_idx].weights.end());
          chunks[b_idx].chars.erase(chunks[b_idx].chars.begin() + mid, chunks[b_idx].chars.end());
          chunks[b_idx].weights.erase(chunks[b_idx].weights.begin() + mid, chunks[b_idx].weights.end());
          chunks.insert(chunks.begin() + b_idx + 1, std::move(new_chunk));
        }
        return;
      }
      cur += chunks[b_idx].size();
    }
  }

  char access(size_t i) const {
    assert(i < total_size);
    size_t cur = 0;
    for (const auto& chk : chunks) {
      if (i < cur + chk.size()) return chk.chars[i - cur];
      cur += chk.size();
    }
    return 0;
  }

  w_t weight_at(size_t i) const {
    assert(i < total_size);
    size_t cur = 0;
    for (const auto& chk : chunks) {
      if (i < cur + chk.size()) return chk.weights[i - cur];
      cur += chk.size();
    }
    return w_t{};
  }

  void set_weight(size_t i, w_t w) {
    assert(i < total_size);
    size_t cur = 0;
    for (auto& chk : chunks) {
      if (i < cur + chk.size()) {
        chk.weights[i - cur] = w;
        return;
      }
      cur += chk.size();
    }
  }

  void range_set_weight(size_t l, size_t r, w_t w) {
    if (l > r || r >= total_size) return;
    size_t cur = 0;
    for (auto& chk : chunks) {
      size_t blk_len = chk.size();
      if (cur + blk_len > l && cur <= r) {
        size_t start = (l > cur) ? (l - cur) : 0;
        size_t end = (r < cur + blk_len - 1) ? (r - cur) : (blk_len - 1);
        for (size_t k = start; k <= end; ++k) {
          chk.weights[k] = w;
        }
      }
      cur += blk_len;
    }
  }

  void remove(size_t i) {
    assert(i < total_size);
    size_t cur = 0;
    for (size_t b_idx = 0; b_idx < chunks.size(); ++b_idx) {
      if (i < cur + chunks[b_idx].size()) {
        size_t offset = i - cur;
        chunks[b_idx].chars.erase(chunks[b_idx].chars.begin() + offset);
        chunks[b_idx].weights.erase(chunks[b_idx].weights.begin() + offset);
        total_size--;
        if (chunks[b_idx].empty() && chunks.size() > 1) {
          chunks.erase(chunks.begin() + b_idx);
        }
        return;
      }
      cur += chunks[b_idx].size();
    }
  }

  size_t size_in_bits() const {
    size_t bits = sizeof(*this) * 8 + chunks.capacity() * sizeof(Chunk) * 8;
    for (const auto& chk : chunks) bits += chk.size_in_bits() - sizeof(Chunk) * 8;
    return bits;
  }
};

template <typename monoid,
          uint32_t B_leaf = 256,
          uint32_t B_fan_out = 16>
class weighted_wt_string {
public:
  using w_t = typename monoid::value_type;
  using Leaf = packed_bit_vector<monoid>;
  using Tree = w_spsi<Leaf, B_leaf, B_fan_out>;
  using BV = weighted_bitvector<Tree, monoid>;

private:
  struct Node {
    uint16_t low;
    uint16_t high;
    BV bv;
    Node* parent = nullptr;
    Node* left = nullptr;
    Node* right = nullptr;
    weighted_wt_string* host = nullptr;

    Node(uint16_t l, uint16_t h, Node* p = nullptr, weighted_wt_string* hst = nullptr)
        : low(l), high(h), parent(p), host(hst) {
      bv.set_resolver(this, [](void* ctx, uint64_t idx) -> w_t {
        auto* self = static_cast<Node*>(ctx);
        return self->host->get_weight_at(self, idx);
      });
    }

    size_t size_in_bits() const {
      size_t bits = sizeof(*this) * 8;
      bits += bv.size_in_bits() - sizeof(bv) * 8;
      if (left) bits += left->size_in_bits();
      if (right) bits += right->size_in_bits();
      return bits;
    }

    ~Node() {
      delete left;
      delete right;
    }
  };

  Node* root_ = nullptr;
  base_leaf_storage<w_t> base_;
  static constexpr uint16_t SIGMA = 256;

  w_t get_weight_at(const Node* node, size_t idx) const {
    const Node* cur = node;
    size_t cur_i = idx;
    while (cur->parent != nullptr) {
      const Node* p = cur->parent;
      if (cur == p->left) {
        cur_i = p->bv.select(0, cur_i + 1);
      } else {
        cur_i = p->bv.select(1, cur_i + 1);
      }
      cur = p;
    }
    return base_.weight_at(cur_i);
  }

  void update_host(Node* n) {
    if (!n) return;
    n->host = this;
    update_host(n->left);
    update_host(n->right);
  }

public:
  weighted_wt_string() {
    root_ = new Node(0, SIGMA - 1, nullptr, this);
  }

  ~weighted_wt_string() {
    delete root_;
  }

  weighted_wt_string(const weighted_wt_string&) = delete;
  weighted_wt_string& operator=(const weighted_wt_string&) = delete;

  weighted_wt_string(weighted_wt_string&& o) noexcept
      : root_(o.root_), base_(std::move(o.base_)) {
    o.root_ = nullptr;
    update_host(root_);
  }

  weighted_wt_string& operator=(weighted_wt_string&& o) noexcept {
    if (this != &o) {
      delete root_;

      root_ = o.root_;
      base_ = std::move(o.base_);
      o.root_ = nullptr;

      update_host(root_);
    }

    return *this;
  }

  size_t size() const {
    return root_ ? root_->bv.size() : 0;
  }

  bool empty() const {
    return size() == 0;
  }

  char access(size_t i) const {
    assert(i < size());
    return base_.access(i);
  }

  w_t weight_at(size_t i) const {
    assert(i < size());
    return base_.weight_at(i);
  }

  void push_back(char c, w_t w) {
    base_.push_back(c, w);
    auto uc = static_cast<unsigned char>(c);
    Node* cur = root_;

    while (true) {
      uint16_t mid = cur->low + (cur->high - cur->low) / 2;
      bool b = (cur->low < cur->high) ? (uc > mid) : false;

      cur->bv.push_back(b, w);
      if (cur->low == cur->high) break;

      if (!b) {
        if (!cur->left) cur->left = new Node(cur->low, mid, cur, this);
        cur = cur->left;
      } else {
        if (!cur->right) cur->right = new Node(mid + 1, cur->high, cur, this);
        cur = cur->right;
      }
    }
  }

  void insert(size_t i, char c, w_t w) {
    assert(i <= size());
    base_.insert(i, c, w);

    auto uc = static_cast<unsigned char>(c);
    Node* cur = root_;

    size_t cur_i = i;

    while (true) {
      uint16_t mid = cur->low + (cur->high - cur->low) / 2;
      bool b = (cur->low < cur->high) ? (uc > mid) : false;

      cur->bv.insert(cur_i, b, w);

      if (cur->low == cur->high) break;

      cur_i = cur->bv.rank(b, cur_i) - 1;

      if (!b) {
        if (!cur->left) cur->left = new Node(cur->low, mid, cur, this);
        cur = cur->left;
      } else {
        if (!cur->right) cur->right = new Node(mid + 1, cur->high, cur, this);
        cur = cur->right;
      }
    }
  }

  void remove(size_t i) {
    assert(i < size());
    base_.remove(i);

    Node* cur = root_;
    size_t cur_i = i;

    while (true) {
      if (cur->low == cur->high) {
        cur->bv.remove(cur_i);
        break;
      }

      bool b = cur->bv.access(cur_i);
      size_t next_i = cur->bv.rank(b, cur_i) - 1;

      cur->bv.remove(cur_i);

      if (!b) {
        assert(cur->left);

        cur = cur->left;
      } else {
        assert(cur->right);

        cur = cur->right;
      }

      cur_i = next_i;
    }
  }

  void set_weight(size_t i, w_t w) {
    assert(i < size());
    base_.set_weight(i, w);

    Node* cur = root_;
    size_t cur_i = i;

    while (true) {
      cur->bv.set_weight(cur_i, w);

      if (cur->low == cur->high) break;

      bool b = cur->bv.access(cur_i);
      cur_i = cur->bv.rank(b, cur_i) - 1;

      if (!b) {
        assert(cur->left);

        cur = cur->left;
      } else {
        assert(cur->right);

        cur = cur->right;
      }
    }
  }

  void range_set_weight(size_t l, size_t r, w_t w) {
    if (l > r || r >= size()) return;
    base_.range_set_weight(l, r, w);
    range_set_weight_impl(root_, l, r, w);
  }

  w_t range_aggregate(size_t l, size_t r, char a, char b) const {
    unsigned char ua = static_cast<unsigned char>(a);
    unsigned char ub = static_cast<unsigned char>(b);

    if (ua > ub || l > r || r >= size() || empty()) {
      return monoid::id;
    }

    return range_aggregate_impl(root_, l, r, ua, ub);
  }

  w_t range_aggregate(size_t l, size_t r) const {
    return range_aggregate(l, r, 0, 255);
  }

  size_t size_in_bits() const {
    size_t bits = sizeof(*this) * 8;
    bits += base_.size_in_bits();
    if (root_) {
      bits += root_->size_in_bits();
    }
    return bits;
  }

  size_t size_in_bytes() const {
    return (size_in_bits() + 7) / 8;
  }

private:
  void range_set_weight_impl(Node* node, size_t l, size_t r, w_t w) {
    if (!node || l > r || node->bv.empty()) return;
    node->bv.range_set_weight(l, r, w);

    if (node->low == node->high) return;

    size_t l0 = (l == 0 ? 0 : node->bv.rank(0, l - 1));
    size_t r0 = node->bv.rank(0, r) - 1;

    if (l0 <= r0 && node->left) {
      range_set_weight_impl(node->left, l0, r0, w);
    }

    size_t l1 = (l == 0 ? 0 : node->bv.rank(1, l - 1));
    size_t r1 = node->bv.rank(1, r) - 1;

    if (l1 <= r1 && node->right) {
      range_set_weight_impl(node->right, l1, r1, w);
    }
  }

  w_t range_aggregate_impl(Node* node, size_t l, size_t r, unsigned char ua, unsigned char ub) const {
    if (!node || l > r || node->bv.empty()) return monoid::id;
    if (node->high < ua || node->low > ub) return monoid::id;

    if (ua <= node->low && node->high <= ub) {
      return node->bv.range_aggregate(l, r);
    }

    uint16_t mid = node->low + (node->high - node->low) / 2;

    size_t l0 = (l == 0 ? 0 : node->bv.rank(0, l - 1));
    size_t r0 = node->bv.rank(0, r) - 1;

    size_t l1 = (l == 0 ? 0 : node->bv.rank(1, l - 1));
    size_t r1 = node->bv.rank(1, r) - 1;

    w_t left_res = monoid::id;
    w_t right_res = monoid::id;

    if (l0 <= r0 && node->left && mid >= ua) {
      left_res = range_aggregate_impl(node->left, l0, r0, ua, ub);
    }

    if (l1 <= r1 && node->right && mid + 1 <= ub) {
      right_res = range_aggregate_impl(node->right, l1, r1, ua, ub);
    }

    return monoid::combine(left_res, right_res);
  }
};

#endif // WEIGHTED_WT_STRING_H