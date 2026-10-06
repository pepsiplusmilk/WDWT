//
// Created by Alex Vladimirov on 06.10.2026.
//

#ifndef TREAP
#define TREAP

#ifndef PRUNED_TREAP_HPP
#define PRUNED_TREAP_HPP

#include <cstdint>
#include <cstddef>
#include <cassert>
#include <algorithm>

template <typename monoid>
class pruned_treap {
public:
    using w_t = typename monoid::value_type;

    struct Node {
        char c;
        w_t w;
        uint32_t priority;
        size_t size = 1;
        w_t subtree_weight;
        uint8_t c_min, c_max;
        Node* left = nullptr;
        Node* right = nullptr;

        Node(char ch, w_t wt, uint32_t prio)
            : c(ch), w(wt), priority(prio), size(1), subtree_weight(wt),
              c_min(static_cast<uint8_t>(ch)), c_max(static_cast<uint8_t>(ch)) {}

        ~Node() { delete left; delete right; }

        size_t size_in_bits() const {
            size_t bits = sizeof(*this) * 8;
            if (left) bits += left->size_in_bits();
            if (right) bits += right->size_in_bits();
            return bits;
        }
    };

private:
    Node* root = nullptr;
    uint32_t rng_state = 13371337;

    uint32_t next_priority() {
        rng_state ^= rng_state << 13; rng_state ^= rng_state >> 17; rng_state ^= rng_state << 5;
        return rng_state;
    }

    static size_t get_size(Node* v) { return v ? v->size : 0; }

    static void update(Node* v) {
        if (!v) return;
        v->size = 1 + get_size(v->left) + get_size(v->right);
        uint8_t uc = static_cast<uint8_t>(v->c);
        v->c_min = v->c_max = uc;

        w_t l_agg = v->left ? v->left->subtree_weight : monoid::id;
        w_t r_agg = v->right ? v->right->subtree_weight : monoid::id;
        v->subtree_weight = monoid::combine(monoid::combine(l_agg, v->w), r_agg);

        if (v->left) {
            v->c_min = std::min(v->c_min, v->left->c_min);
            v->c_max = std::max(v->c_max, v->left->c_max);
        }
        if (v->right) {
            v->c_min = std::min(v->c_min, v->right->c_min);
            v->c_max = std::max(v->c_max, v->right->c_max);
        }
    }

    static void split(Node* t, size_t k, Node*& l, Node*& r) {
        if (!t) { l = r = nullptr; return; }
        size_t left_sz = get_size(t->left);
        if (k <= left_sz) {
            split(t->left, k, l, t->left);
            r = t;
        } else {
            split(t->right, k - left_sz - 1, t->right, r);
            l = t;
        }
        update(t);
    }

    static void merge(Node*& t, Node* l, Node* r) {
        if (!l || !r) { t = l ? l : r; return; }
        if (l->priority > r->priority) {
            merge(l->right, l->right, r);
            t = l;
        } else {
            merge(r->left, l, r->left);
            t = r;
        }
        update(t);
    }

    static void set_weight_impl(Node* cur, size_t i, w_t w) {
        if (!cur) return;
        size_t l_sz = get_size(cur->left);
        if (i < l_sz) set_weight_impl(cur->left, i, w);
        else if (i == l_sz) cur->w = w;
        else set_weight_impl(cur->right, i - l_sz - 1, w);
        update(cur);
    }

    static w_t query_pruned(Node* v, uint8_t ua, uint8_t ub) {
        if (!v || v->c_max < ua || v->c_min > ub) return monoid::id;
        if (ua <= v->c_min && v->c_max <= ub) return v->subtree_weight;

        w_t res = monoid::id;
        uint8_t uc = static_cast<uint8_t>(v->c);
        if (ua <= uc && uc <= ub) res = v->w;
        if (v->left) res = monoid::combine(query_pruned(v->left, ua, ub), res);
        if (v->right) res = monoid::combine(res, query_pruned(v->right, ua, ub));
        return res;
    }

public:
    pruned_treap() = default;
    ~pruned_treap() { delete root; }

    size_t size() const { return get_size(root); }
    bool empty() const { return size() == 0; }

    void insert(size_t i, char c, w_t w) {
        Node *l, *r;
        split(root, i, l, r);
        Node* m = new Node(c, w, next_priority());
        merge(l, l, m); merge(root, l, r);
    }

    void push_back(char c, w_t w) { insert(size(), c, w); }

    void remove(size_t i) {
        Node *l, *m, *r;
        split(root, i, l, m); split(m, 1, m, r);
        delete m;
        merge(root, l, r);
    }

    char access(size_t i) const {
        Node* cur = root;
        while (cur) {
            size_t l_sz = get_size(cur->left);
            if (i < l_sz) cur = cur->left;
            else if (i == l_sz) return cur->c;
            else { i -= l_sz + 1; cur = cur->right; }
        }
        return 0;
    }

    w_t weight_at(size_t i) const {
        Node* cur = root;
        while (cur) {
            size_t l_sz = get_size(cur->left);
            if (i < l_sz) cur = cur->left;
            else if (i == l_sz) return cur->w;
            else { i -= l_sz + 1; cur = cur->right; }
        }
        return monoid::id;
    }

    void set_weight(size_t i, w_t w) { set_weight_impl(root, i, w); }

    void range_set_weight(size_t l, size_t r, w_t w) {
        if (l > r || r >= size()) return;
        for (size_t i = l; i <= r; ++i) set_weight(i, w);
    }

    w_t range_aggregate(size_t l, size_t r, char a, char b) {
        unsigned char ua = static_cast<unsigned char>(a), ub = static_cast<unsigned char>(b);
        if (ua > ub || l > r || r >= size() || empty()) return monoid::id;

        Node *t1, *t2, *t3;
        split(root, l, t1, t2); split(t2, r - l + 1, t2, t3);
        w_t ans = query_pruned(t2, ua, ub);
        merge(t2, t2, t3); merge(root, t1, t2);
        return ans;
    }

    w_t range_aggregate(size_t l, size_t r) { return range_aggregate(l, r, 0, 255); }

    size_t size_in_bits() const { return sizeof(*this) * 8 + (root ? root->size_in_bits() : 0); }
    size_t size_in_bytes() const { return (size_in_bits() + 7) / 8; }
};
#endif // PRUNED_TREAP_HPP

#endif  // TREAP