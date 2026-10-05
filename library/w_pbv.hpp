//
// Created by Alex Vladimirov on 05.10.2026.
//

#ifndef W_PBV_H
#define W_PBV_H

#include<vector>
#include <cassert>
#include <cstdint>

template<typename monoid>
class packed_bit_vector {
public:
  using w_t = typename monoid::value_type;

  explicit packed_bit_vector(size_t size = 0, bool neutral_bit = 0, w_t neutral_w = monoid::id) {
    if (size > 0) {
      size_ = size;

      size_t size_in_words = (size + 63) / 64;
      words.assign(size_in_words, neutral_bit ? ~static_cast<uint64_t>(0) : static_cast<uint64_t>(0));

      if (neutral_bit && (size % 64 != 0)) {
        words.back() &= (uint64_t(1) << (size % 64)) - 1;
      }

      weights.assign(size, neutral_w);

      rebuild_aggregate();
    }
  }

  packed_bit_vector(std::vector<uint64_t>&& w, std::vector<w_t>&& wt, uint64_t sz)
      : words(std::move(w)), weights(std::move(wt)), size_(sz) {
    rebuild_aggregate();
  }

  uint64_t size() const {
    return size_;
  }

  uint64_t psum() const {
    return psum_;
  }

  w_t weight() const {
    return weight_aggregate_;
  }

  bool at(uint64_t i) const {
    assert(i < size_);

    return (words[i / 64] >> (i % 64)) & 1ull;
  }

  const w_t& weight_at(uint64_t i) const {
    assert(i < size_);

    return weights[i];
  }

  uint64_t psum(uint64_t i) const {
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
        uint64_t cur_size = std::min<uint64_t>(64, size_ - w * 64);
        bits = cur_size - bits;
      }

      if (s + bits >= x) break;
      s += bits;
      w++;
    }

    uint64_t word = words[w];
    uint64_t bit_pos = 0;
    while (s < x) {
      if (((word >> bit_pos) & 1ull) == b) {
        s++;
      }
      if (s == x) break;
      bit_pos++;
    }
    return w * 64 + bit_pos;
  }

  w_t range_aggregate(uint64_t l, uint64_t r) const {
    assert(l <= r && r < size_);

    w_t acc = monoid::id;
    for (uint64_t k = l; k <= r; ++k) {
      acc = monoid::combine(acc, weights[k]);
    }

    return acc;
  }

  void push_back(bool bit, w_t w) {
    if (size_ >= words.size() * 64) {
      words.push_back(0);
    }

    if (bit) {
      words[size_ / 64] |= (1ull << (size_ % 64));
      psum_++;
    }

    weights.push_back(w);
    weight_aggregate_ = monoid::combine(weight_aggregate_, w);
    size_++;
  }

  void insert(uint64_t i, bool bit, w_t w) {
    assert(i <= size_);

    if (i == size_) {
      push_back(bit, w);
      return;
    }

    if (size_ >= words.size() * 64) {
      words.push_back(0);
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

    weights.insert(weights.begin() + i, w);
    size_++;

    rebuild_aggregate();
  }

  void remove(uint64_t i) {
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

    weights.erase(weights.begin() + i);
    size_--;

    if (words.size() * 64 >= size_ + 64 && words.size() > 1) {
      words.pop_back();
    }

    rebuild_aggregate();
  }

  void set(uint64_t i, bool bit, w_t w) {
    assert(i < size_);

    bool old_bit = at(i);

    if (old_bit != bit) {
      if (bit) {
        words[i / 64] |= (1ull << (i % 64));
      } else {
        words[i / 64] &= ~(1ull << (i % 64));
      }
    }

    weights[i] = w;
    rebuild_aggregate();
  }

  packed_bit_vector* split() {
    uint64_t tot_words = (size_ + 63) / 64;

    assert(tot_words >= 2);

    uint64_t nr_left_words = tot_words / 2;
    uint64_t nr_left_ints = nr_left_words * 64;
    uint64_t nr_right_ints = size_ - nr_left_ints;

    std::vector<uint64_t> right_words(words.begin() + nr_left_words, words.end());
    std::vector<w_t> right_weights(weights.begin() + nr_left_ints, weights.end());

    words.resize(nr_left_words);
    weights.resize(nr_left_ints);
    size_ = nr_left_ints;

    auto right = new packed_bit_vector(std::move(right_words), std::move(right_weights), nr_right_ints);

    rebuild_aggregate();

    return right;
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
    for (const auto& w : weights) {
      weight_aggregate_ = monoid::combine(weight_aggregate_, w);
    }
  }

  std::vector<uint64_t> words;
  std::vector<w_t> weights;

  w_t weight_aggregate_ = monoid::id;

  uint64_t psum_= 0; // aggregate of bits setted to 1 at prefix
  uint64_t size_= 0; // size of leaf
  uint8_t width_= 64; // width of bits that stored as block
  uint8_t int_per_word_= 1;

  //when reallocating, reserve extra_ words of space to accelerate insert
  static const uint8_t extra_ = 2;
};

#endif  // W_PBV_H
