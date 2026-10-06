//
// Created by Alex Vladimirov on 06.10.2026.
//

#ifndef SQRT
#define SQRT

#ifndef SQRT_DECOMPOSITION_HPP
#define SQRT_DECOMPOSITION_HPP

#include <vector>
#include <array>
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <algorithm>

template <typename monoid, size_t BLOCK_SIZE = 1024>
class sqrt_decomposition {
public:
    using w_t = typename monoid::value_type;

    struct Block {
        std::vector<char> chars;
        std::vector<w_t> weights;
        std::array<w_t, 256> block_sum;

        Block() {
            chars.reserve(BLOCK_SIZE * 2);
            weights.reserve(BLOCK_SIZE * 2);
            block_sum.fill(monoid::id);
        }

        void rebuild() {
            block_sum.fill(monoid::id);
            for (size_t i = 0; i < chars.size(); ++i) {
                uint8_t c = static_cast<uint8_t>(chars[i]);
                block_sum[c] = monoid::combine(block_sum[c], weights[i]);
            }
        }

        size_t size() const { return chars.size(); }
        bool empty() const { return chars.empty(); }

        size_t size_in_bits() const {
            size_t bits = sizeof(*this) * 8;
            bits += chars.capacity() * sizeof(char) * 8;
            bits += weights.capacity() * sizeof(w_t) * 8;
            return bits;
        }
    };

private:
    std::vector<Block> blocks;
    size_t total_size = 0;

public:
    sqrt_decomposition() { blocks.emplace_back(); }
    size_t size() const { return total_size; }
    bool empty() const { return total_size == 0; }

    void push_back(char c, w_t w) { insert(total_size, c, w); }

    void insert(size_t i, char c, w_t w) {
        assert(i <= total_size);
        size_t cur = 0, b_idx = 0;
        for (; b_idx < blocks.size(); ++b_idx) {
            if (i <= cur + blocks[b_idx].size()) {
                if (b_idx + 1 == blocks.size() || i < cur + blocks[b_idx].size()) break;
            }
            cur += blocks[b_idx].size();
        }
        if (b_idx == blocks.size()) {
            b_idx = blocks.size() - 1;
            cur = total_size - blocks.back().size();
        }

        size_t offset = i - cur;
        Block& blk = blocks[b_idx];
        blk.chars.insert(blk.chars.begin() + offset, c);
        blk.weights.insert(blk.weights.begin() + offset, w);
        blk.rebuild();
        total_size++;

        if (blk.size() >= 2 * BLOCK_SIZE) {
            Block new_blk;
            size_t mid = blk.size() / 2;
            new_blk.chars.assign(blk.chars.begin() + mid, blk.chars.end());
            new_blk.weights.assign(blk.weights.begin() + mid, blk.weights.end());
            blk.chars.erase(blk.chars.begin() + mid, blk.chars.end());
            blk.weights.erase(blk.weights.begin() + mid, blk.weights.end());
            blk.rebuild();
            new_blk.rebuild();
            blocks.insert(blocks.begin() + b_idx + 1, std::move(new_blk));
        }
    }

    char access(size_t i) const {
        size_t cur = 0;
        for (const auto& blk : blocks) {
            if (i < cur + blk.size()) return blk.chars[i - cur];
            cur += blk.size();
        }
        return 0;
    }

    w_t weight_at(size_t i) const {
        size_t cur = 0;
        for (const auto& blk : blocks) {
            if (i < cur + blk.size()) return blk.weights[i - cur];
            cur += blk.size();
        }
        return monoid::id;
    }

    void remove(size_t i) {
        size_t cur = 0;
        for (size_t b_idx = 0; b_idx < blocks.size(); ++b_idx) {
            if (i < cur + blocks[b_idx].size()) {
                size_t offset = i - cur;
                blocks[b_idx].chars.erase(blocks[b_idx].chars.begin() + offset);
                blocks[b_idx].weights.erase(blocks[b_idx].weights.begin() + offset);
                blocks[b_idx].rebuild();
                total_size--;
                if (blocks[b_idx].empty() && blocks.size() > 1) {
                    blocks.erase(blocks.begin() + b_idx);
                }
                return;
            }
            cur += blocks[b_idx].size();
        }
    }

    void set_weight(size_t i, w_t w) {
        size_t cur = 0;
        for (auto& blk : blocks) {
            if (i < cur + blk.size()) {
                blk.weights[i - cur] = w;
                blk.rebuild();
                return;
            }
            cur += blk.size();
        }
    }

    void range_set_weight(size_t l, size_t r, w_t w) {
        if (l > r || r >= total_size) return;
        size_t cur = 0;
        for (auto& blk : blocks) {
            size_t blk_len = blk.size();
            if (cur + blk_len > l && cur <= r) {
                size_t start = (l > cur) ? (l - cur) : 0;
                size_t end = (r < cur + blk_len - 1) ? (r - cur) : (blk_len - 1);
                for (size_t k = start; k <= end; ++k) blk.weights[k] = w;
                blk.rebuild();
            }
            cur += blk_len;
        }
    }

    w_t range_aggregate(size_t l, size_t r, char a, char b) const {
        unsigned char ua = static_cast<unsigned char>(a), ub = static_cast<unsigned char>(b);
        if (ua > ub || l > r || r >= total_size || empty()) return monoid::id;

        w_t res = monoid::id;
        size_t cur = 0;
        for (const auto& blk : blocks) {
            size_t blk_len = blk.size();
            if (cur + blk_len > l && cur <= r) {
                size_t start = (l > cur) ? (l - cur) : 0;
                size_t end = (r < cur + blk_len - 1) ? (r - cur) : (blk_len - 1);
                if (start == 0 && end == blk_len - 1) {
                    for (uint16_t c = ua; c <= ub; ++c) res = monoid::combine(res, blk.block_sum[c]);
                } else {
                    for (size_t k = start; k <= end; ++k) {
                        uint8_t c = static_cast<uint8_t>(blk.chars[k]);
                        if (c >= ua && c <= ub) res = monoid::combine(res, blk.weights[k]);
                    }
                }
            }
            cur += blk_len;
        }
        return res;
    }

    w_t range_aggregate(size_t l, size_t r) const { return range_aggregate(l, r, 0, 255); }

    size_t size_in_bits() const {
        size_t bits = sizeof(*this) * 8 + blocks.capacity() * sizeof(Block) * 8;
        for (const auto& blk : blocks) bits += blk.size_in_bits() - sizeof(Block) * 8;
        return bits;
    }

    size_t size_in_bytes() const { return (size_in_bits() + 7) / 8; }
};
#endif // SQRT_DECOMPOSITION_HPP

#endif  // SQRT