//
// Created by Alex Vladimirov on 05.10.2026.
//

#ifndef W_PBV_H
#define W_PBV_H

#include <vector>
#include <cassert>
#include <cstdint>
#include <algorithm>

#if defined(__BMI2__)
#include <immintrin.h>
#endif

template <typename w_t>
struct weight_resolver {
  void* ctx = nullptr;
  using fn_t = w_t (*)(void*, uint64_t);
  fn_t fn = nullptr;

  inline w_t get(uint64_t idx) const {
    return fn ? fn(ctx, idx) : w_t{};
  }
  explicit operator bool() const { return fn != nullptr; }
};

// vectorized select
static inline uint64_t select64(uint64_t word, uint64_t k) {
#if defined(__BMI2__)
  return __builtin_ctzll(_pdep_u64(1ull << (k - 1), word));
#else
  uint64_t pos = 0;

  uint64_t pop32 = __builtin_popcountll(word & 0xFFFFFFFFull);
  if (k > pop32) { k -= pop32; pos += 32; word >>= 32; }

  uint64_t pop16 = __builtin_popcountll(word & 0xFFFFull);
  if (k > pop16) { k -= pop16; pos += 16; word >>= 16; }

  uint64_t pop8 = __builtin_popcountll(word & 0xFFull);
  if (k > pop8) { k -= pop8; pos += 8; word >>= 8; }

  while (k > 0) {
    if (word & 1ull) {
      k--;

      if (k == 0) {
        return pos;
      }
    }

    word >>= 1;
    pos++;
  }

  return pos;
#endif
}

template<typename monoid>
class packed_bit_vector {
public:
  using w_t = typename monoid::value_type;
  using monoid_type = monoid;
  using resolver_t = weight_resolver<w_t>;

  explicit packed_bit_vector(size_t size = 0, bool neutral_bit = 0, w_t neutral_w = monoid::id) {
    if (size > 0) {
      size_ = size;

      size_t size_in_words = (size + 63) / 64;
      words.assign(size_in_words, neutral_bit ? ~static_cast<uint64_t>(0) : static_cast<uint64_t>(0));

      if (neutral_bit && (size % 64 != 0)) {
        words.back() &= (uint64_t(1) << (size % 64)) - 1;
      }

      block_weights.assign(size_in_words, neutral_w);
      rebuild_aggregate();
    }
  }

  packed_bit_vector(std::vector<uint64_t>&& w, std::vector<w_t>&& bw, uint64_t sz)
      : words(std::move(w)), block_weights(std::move(bw)), size_(sz) {
    rebuild_aggregate();
  }

  uint64_t size() const {
    return size_;
  }

  bool empty() const {
    return size_ == 0;
  }

  uint64_t psum() const {
    return psum_;
  }

  inline bool at(uint64_t i) const {
    assert(i < size_);

    return (words[i / 64] >> (i % 64)) & 1ull;
  }

  inline uint64_t psum(uint64_t i) const {
    assert(i < size_);

    uint64_t s = 0;
    size_t word_idx = i / 64;

    for (size_t j = 0; j < word_idx; ++j) {
      s += __builtin_popcountll(words[j]);
    }

    uint8_t bit_offset = i % 64;
    uint64_t mask = (bit_offset == 63) ? ~static_cast<uint64_t>(0) : ((static_cast<uint64_t>(1) << (bit_offset + 1)) - 1);
    s += __builtin_popcountll(words[word_idx] & mask);

    return s;
  }

  uint64_t search(uint64_t x, bool b) const {
    uint64_t count = b ? psum_ : size_ - psum_;
    assert(size_ > 0 && x <= count);

    if (!x) {
      return 0;
    }

    uint64_t s = 0;
    size_t w = 0;
    while (w < words.size()) {
      uint64_t bits = __builtin_popcountll(words[w]);

      if (!b) {
        uint64_t cur_size = std::min<uint64_t>(64, size_ > w * 64 ? size_ - w * 64 : 0);
        bits = cur_size - bits;
      }

      if (s + bits >= x) break;
      s += bits;
      w++;
    }

    uint64_t word = (w < words.size() ? words[w] : 0);
    uint64_t needed = x - s;
    uint64_t bit_pos = b ? select64(word, needed) : select64(~word, needed);
    return w * 64 + bit_pos;
  }

  w_t range_aggregate(uint64_t l, uint64_t r, resolver_t res, uint64_t leaf_offset) const {
    assert(l <= r && r < size_);
    size_t start_w = l / 64;
    size_t end_w = r / 64;

    if (start_w == end_w) {
      if (l == start_w * 64 && (r == start_w * 64 + 63 || r == size_ - 1)) {
        return block_weights[start_w];
      }
      w_t acc = monoid::id;
      for (uint64_t k = l; k <= r; ++k) {
        acc = monoid::combine(acc, res.get(leaf_offset + k));
      }
      return acc;
    }

    w_t acc = monoid::id;
    if (l == start_w * 64) {
      acc = monoid::combine(acc, block_weights[start_w]);
    } else {
      for (uint64_t k = l; k < (start_w + 1) * 64; ++k) {
        acc = monoid::combine(acc, res.get(leaf_offset + k));
      }
    }

    for (size_t w = start_w + 1; w < end_w; ++w) {
      acc = monoid::combine(acc, block_weights[w]);
    }

    uint64_t end_w_last = std::min(size_ - 1, static_cast<uint64_t>(end_w * 64 + 63));
    if (r == end_w_last) {
      acc = monoid::combine(acc, block_weights[end_w]);
    } else {
      for (uint64_t k = end_w * 64; k <= r; ++k) {
        acc = monoid::combine(acc, res.get(leaf_offset + k));
      }
    }

    return acc;
  }

  void set_weight(uint64_t i, w_t /*w*/, resolver_t res, uint64_t leaf_offset) {
    assert(i < size_);

    size_t w_idx = i / 64;
    uint64_t start_k = w_idx * 64;
    uint64_t end_k = std::min(size_, static_cast<uint64_t>((w_idx + 1) * 64));

    w_t acc = monoid::id;
    for (uint64_t k = start_k; k < end_k; ++k) {
      acc = monoid::combine(acc, res.get(leaf_offset + k));
    }
    block_weights[w_idx] = acc;

    rebuild_aggregate();
  }

  void range_set_weight(uint64_t l, uint64_t r, w_t /*w*/, resolver_t res, uint64_t leaf_offset) {
    assert(l <= r && r < size_);
    size_t start_w = l / 64;
    size_t end_w = r / 64;

    for (size_t w_idx = start_w; w_idx <= end_w; ++w_idx) {
      uint64_t start_k = w_idx * 64;
      uint64_t end_k = std::min(size_, static_cast<uint64_t>((w_idx + 1) * 64));
      w_t acc = monoid::id;
      for (uint64_t k = start_k; k < end_k; ++k) {
        acc = monoid::combine(acc, res.get(leaf_offset + k));
      }
      block_weights[w_idx] = acc;
    }

    rebuild_aggregate();
  }

  w_t total_weight() const {
    return weight_aggregate_;
  }

  void push_back(bool bit, w_t w) {
    if (size_ >= words.size() * 64) {
      words.push_back(0);
      block_weights.push_back(monoid::id);
    }

    if (bit) {
      words[size_ / 64] |= (1ull << (size_ % 64));
      psum_++;
    }

    size_t w_idx = size_ / 64;
    block_weights[w_idx] = monoid::combine(block_weights[w_idx], w);
    weight_aggregate_ = monoid::combine(weight_aggregate_, w);
    size_++;
  }

  void insert(uint64_t i, bool bit, w_t w, resolver_t res, uint64_t leaf_offset) {
    assert(i <= size_);

    if (i == size_) {
      push_back(bit, w);
      return;
    }

    if (size_ >= words.size() * 64) {
      words.push_back(0);
      block_weights.push_back(monoid::id);
    }

    size_t start_w = i / 64;
    size_t bit_pos = i % 64;

    uint64_t carry = bit ? 1ull : 0ull;

    uint64_t mask = (bit_pos == 0) ? 0ull : ((1ull << bit_pos) - 1ull);

    uint64_t low = words[start_w] & mask;
    uint64_t high = words[start_w] & ~mask;

    uint64_t next_carry = (words[start_w] >> 63) & 1ull;

    words[start_w] = low | (carry << bit_pos) | (high << 1);
    carry = next_carry;

    for (size_t w_idx = start_w + 1; w_idx < words.size(); ++w_idx) {
      next_carry = (words[w_idx] >> 63) & 1ull;
      words[w_idx] = (words[w_idx] << 1) | carry;

      carry = next_carry;
    }

    size_++;

    for (size_t w_idx = start_w; w_idx < words.size(); ++w_idx) {
      uint64_t start_k = w_idx * 64;
      uint64_t end_k = std::min(size_, static_cast<uint64_t>((w_idx + 1) * 64));
      w_t acc = monoid::id;
      for (uint64_t k = start_k; k < end_k; ++k) {
        acc = monoid::combine(acc, res.get(leaf_offset + k));
      }
      block_weights[w_idx] = acc;
    }

    rebuild_aggregate();
  }

  void remove(uint64_t i, resolver_t res, uint64_t leaf_offset) {
    assert(i < size_);

    size_t start_w = i / 64;
    size_t bit_pos = i % 64;

    uint64_t mask = (bit_pos == 0) ? 0ull : ((1ull << bit_pos) - 1ull);

    uint64_t low = words[start_w] & mask;
    uint64_t high = (words[start_w] >> 1) & ~mask;

    words[start_w] = low | high;

    for (size_t w_idx = start_w + 1; w_idx < words.size(); ++w_idx) {
      words[w_idx - 1] |= (words[w_idx] & 1ull) << 63;
      words[w_idx] >>= 1;
    }

    size_--;

    if (words.size() * 64 >= size_ + 64 && words.size() > 1) {
      words.pop_back();
      block_weights.pop_back();
    }

    for (size_t w_idx = start_w; w_idx < words.size(); ++w_idx) {
      uint64_t start_k = w_idx * 64;
      uint64_t end_k = std::min(size_, static_cast<uint64_t>((w_idx + 1) * 64));
      w_t acc = monoid::id;
      for (uint64_t k = start_k; k < end_k; ++k) {
        acc = monoid::combine(acc, res.get(leaf_offset + k));
      }
      block_weights[w_idx] = acc;
    }

    rebuild_aggregate();
  }

  packed_bit_vector* split() {
    uint64_t tot_words = (size_ + 63) / 64;
    uint64_t nr_left_words = tot_words / 2;
    uint64_t nr_left_ints = nr_left_words * 64;
    uint64_t nr_right_ints = size_ - nr_left_ints;

    if (tot_words >= 2) {
      std::vector<uint64_t> right_words(words.begin() + nr_left_words, words.end());
      std::vector<w_t> right_bw(block_weights.begin() + nr_left_words, block_weights.end());

      words.resize(nr_left_words);
      block_weights.resize(nr_left_words);
      size_ = nr_left_ints;

      auto right = new packed_bit_vector(std::move(right_words), std::move(right_bw), nr_right_ints);
      rebuild_aggregate();
      return right;
    } else {
      nr_left_ints = size_ / 2;
      nr_right_ints = size_ - nr_left_ints;

      auto right = new packed_bit_vector();
      for (uint64_t k = nr_left_ints; k < size_; ++k) {
        right->push_back(at(k), block_weights.empty() ? monoid::id : block_weights[0]);
      }

      words.resize(1);
      words[0] &= (nr_left_ints == 64 ? ~0ull : (1ull << nr_left_ints) - 1);
      block_weights.resize(1);
      size_ = nr_left_ints;
      rebuild_aggregate();
      return right;
    }
  }

  size_t size_in_bits() const {
    size_t bits = sizeof(*this) * 8;
    bits += words.capacity() * sizeof(uint64_t) * 8;
    bits += block_weights.capacity() * sizeof(w_t) * 8;
    return bits;
  }

  size_t size_in_bytes() const {
    return (size_in_bits() + 7) / 8;
  }

  virtual ~packed_bit_vector() = default;

private:
  void rebuild_aggregate() {
    psum_ = 0;
    for (size_t w = 0; w < words.size(); ++w) {
      uint64_t cur_size = std::min<uint64_t>(64, size_ > w * 64 ? size_ - w * 64 : 0);

      if (cur_size == 64) {
        psum_ += __builtin_popcountll(words[w]);
      } else if (cur_size > 0) {
        uint64_t mask = (1ull << cur_size) - 1;

        psum_ += __builtin_popcountll(words[w] & mask);
      }
    }

    weight_aggregate_ = monoid::id;
    for (const auto& bw : block_weights) {
      weight_aggregate_ = monoid::combine(weight_aggregate_, bw);
    }
  }

  std::vector<uint64_t> words;
  std::vector<w_t> block_weights;

  w_t weight_aggregate_ = monoid::id;
  uint64_t psum_ = 0;
  uint64_t size_ = 0;
};

#endif  // W_PBV_H