//
// Created by Alex Vladimirov on 05.10.2026.
//

#ifndef W_SPSI_NODE_H
#define W_SPSI_NODE_H

#include <array>
#include <cassert>
#include <vector>
#include <algorithm>

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

  ~node() = default;

  uint64_t size() const {
    return nr_children_ ? subtree_sizes[nr_children_ - 1] : 0;
  }

  uint64_t total_ones() const {
    return nr_children_ ? subtree_psums[nr_children_ - 1] : 0;
  }

  w_t total_weight() const {
    return nr_children_ ? subtree_aggs[nr_children_ - 1] : monoid::id;
  }

  bool at(uint64_t i) const {
    assert(i < size());

    uint32_t j = find_child(i);
    uint64_t prev_sz = (j == 0 ? 0 : subtree_sizes[j - 1]);

    if (has_leaves_) {
      return leaves[j]->at(i - prev_sz);
    }

    return children[j]->at(i - prev_sz);
  }

  w_t weight_at(uint64_t i) const {
    assert(i < size());

    uint32_t j = find_child(i);
    uint64_t prev_sz = (j == 0 ? 0 : subtree_sizes[j - 1]);

    if (has_leaves_) {
      return leaves[j]->weight_at(i - prev_sz);
    }

    return children[j]->weight_at(i - prev_sz);
  }

  uint64_t rank_1(uint64_t i) const {
    assert(i < size());

    uint32_t j = find_child(i);

    uint64_t prev_sz = (j == 0 ? 0 : subtree_sizes[j - 1]);
    uint64_t prev_ones = (j == 0 ? 0 : subtree_psums[j - 1]);

    if (has_leaves_) {
      return prev_ones + leaves[j]->psum(i - prev_sz);
    }

    return prev_ones + children[j]->rank_1(i - prev_sz);
  }

  uint64_t rank_0(uint64_t i) const {
    return i + 1 - rank_1(i);
  }

  uint64_t select_1(uint64_t x) const {
    assert(x > 0 && x <= total_ones());

    uint32_t j = find_1(x);

    uint64_t prev_sz = (j == 0 ? 0 : subtree_sizes[j - 1]);
    uint64_t prev_ones = (j == 0 ? 0 : subtree_psums[j - 1]);

    if (has_leaves_) {
      return prev_sz + leaves[j]->search(x - prev_ones, true);
    }

    return prev_sz + children[j]->select_1(x - prev_ones);
  }

  uint64_t select_0(uint64_t x) const {
    assert(x > 0 && x <= (size() - total_ones()));

    uint32_t j = find_0(x);

    uint64_t prev_sz = (j == 0 ? 0 : subtree_sizes[j - 1]);
    uint64_t prev_zeros = prev_sz - (j == 0 ? 0 : subtree_psums[j - 1]);

    if (has_leaves_) {
      return prev_sz + leaves[j]->search(x - prev_zeros, false);
    }

    return prev_sz + children[j]->select_0(x - prev_zeros);
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
      res = monoid::combine(res, mid_weight);
    }

    uint64_t prev_r = subtree_sizes[jr - 1];
    w_t right_res = (has_leaves_ ? leaves[jr]->range_aggregate(0, r - prev_r)
                                 : children[jr]->range_aggregate(0, r - prev_r));

    return monoid::combine(res, right_res);
  }

  void set_weight(uint64_t i, w_t w) {
    assert(i < size());

    uint32_t j = find_child(i);
    uint64_t prev_sz = (j == 0 ? 0 : subtree_sizes[j - 1]);

    if (has_leaves_) {
      leaves[j]->set_weight(i - prev_sz, w);
    } else {
      children[j]->set_weight(i - prev_sz, w);
    }

    rebuild(j);
  }

  void range_set_weight(uint64_t l, uint64_t r, w_t w) {
    assert(l <= r && r < size());

    uint32_t jl = find_child(l);
    uint32_t jr = find_child(r);

    uint64_t prev_l = (jl == 0 ? 0 : subtree_sizes[jl - 1]);

    if (jl == jr) {
      if (has_leaves_) {
        leaves[jl]->range_set_weight(l - prev_l, r - prev_l, w);
      } else {
        children[jl]->range_set_weight(l - prev_l, r - prev_l, w);
      }

      rebuild(jl);
      return;
    }

    uint64_t child_l_sz = (has_leaves_ ? leaves[jl]->size() : children[jl]->size());

    if (has_leaves_) {
      leaves[jl]->range_set_weight(l - prev_l, child_l_sz - 1, w);
    } else {
      children[jl]->range_set_weight(l - prev_l, child_l_sz - 1, w);
    }

    for (uint32_t k = jl + 1; k < jr; ++k) {
      uint64_t sz_k = (has_leaves_ ? leaves[k]->size() : children[k]->size());

      if (sz_k > 0) {
        if (has_leaves_) {
          leaves[k]->range_set_weight(0, sz_k - 1, w);
        } else {
          children[k]->range_set_weight(0, sz_k - 1, w);
        }
      }
    }

    uint64_t prev_r = subtree_sizes[jr - 1];

    if (has_leaves_) {
      leaves[jr]->range_set_weight(0, r - prev_r, w);
    } else {
      children[jr]->range_set_weight(0, r - prev_r, w);
    }

    rebuild(jl);
  }

  node* insert(uint64_t i, bool bit, w_t weight) {
    assert(i <= size());
    node* new_root = nullptr;

    if (is_full()) {
      node* right = split();

      if (i < size()) {
        insert_without_split(i, bit, weight);
      } else {
        right->insert_without_split(i - size(), bit, weight);
      }

      if (is_root()) {
        new_root = new node(std::vector<node*>{this, right});

        this->set_parent(new_root);
        right->set_parent(new_root);
      } else {
        parent_->new_children(rank_, this, right);
      }
    } else {
      insert_without_split(i, bit, weight);
    }

    return new_root;
  }

  node* remove(uint64_t i) {
    assert(i < size());
    node* x = this;

    if (!x->can_loose_child() && !x->is_root()) {
      assert(x->parent_ != nullptr);

      node* y;
      bool y_is_prev;

      if (x->rank_ > 0) {
        y = x->parent_->children[x->rank_ - 1];
        y_is_prev = true;
      } else {
        y = x->parent_->children[x->rank_ + 1];
        y_is_prev = false;
      }

      if (y->can_loose_child()) {
        if (!x->has_leaves_) {
          node* z;
          if (y_is_prev) {
            z = y->children.back();

            y->children.pop_back();
            y->nr_children_--;

            z->set_parent(x);

            x->children.insert(x->children.begin(), z);
            x->nr_children_++;

            i += z->size();
          } else {
            z = y->children.front();

            y->children.erase(y->children.begin());
            y->nr_children_--;

            z->set_parent(x);

            x->children.push_back(z);
            x->nr_children_++;
          }

          uint32_t r = 0;
          for (auto* c : x->children) c->set_rank(r++);
          r = 0;
          for (auto* c : y->children) c->set_rank(r++);

          x->rebuild(0);
          y->rebuild(0);

          x->parent_->rebuild(0);
        } else {
          leaf_type* z;

          if (y_is_prev) {
            z = y->leaves.back();

            y->leaves.pop_back();
            y->nr_children_--;

            x->leaves.insert(x->leaves.begin(), z);
            x->nr_children_++;

            i += z->size();
          } else {
            z = y->leaves.front();

            y->leaves.erase(y->leaves.begin());
            y->nr_children_--;

            x->leaves.push_back(z);
            x->nr_children_++;
          }
          x->rebuild(0);
          y->rebuild(0);

          x->parent_->rebuild(0);
        }
      } else {
        node* prev = y_is_prev ? y : x;
        node* next = y_is_prev ? x : y;

        if (y_is_prev) {
          i += y->size();
        }

        if (!x->has_leaves_) {
          for (auto* c : next->children) {
            c->set_parent(prev);
            prev->children.push_back(c);
          }

          next->children.clear();
        } else {
          for (auto* l : next->leaves) {
            prev->leaves.push_back(l);
          }

          next->leaves.clear();
        }

        prev->nr_children_ += next->nr_children_;
        next->nr_children_ = 0;

        uint32_t next_rank = next->rank_;
        node* P = prev->parent_;

        P->children.erase(P->children.begin() + next_rank);
        P->nr_children_--;

        for (uint32_t r = next_rank; r < P->nr_children_; ++r) {
          P->children[r]->set_rank(r);
        }

        delete next;

        uint32_t r = 0;
        if (!prev->has_leaves_) {
          for (auto* c : prev->children) c->set_rank(r++);
        }

        prev->rebuild(0);
        P->rebuild(0);

        x = prev;
      }
    }

    assert(x->can_loose_child() || x->is_root());

    uint32_t j = x->find_child(i);
    uint64_t prev_sz = (j == 0 ? 0 : x->subtree_sizes[j - 1]);
    uint64_t child_i = i - prev_sz;

    if (x->has_leaves_) {
      leaf_type* lf = x->leaves[j];

      if (!x->leaf_can_loose(lf) && x->leaves.size() > 1) {
        leaf_type* y_leaf = nullptr;
        bool y_is_prev = false;
        if (j > 0) {
          y_leaf = x->leaves[j - 1];
          y_is_prev = true;
        } else {
          y_leaf = x->leaves[j + 1];
          y_is_prev = false;
        }

        if (x->leaf_can_loose(y_leaf)) {
          if (y_is_prev) {
            bool z_bit = y_leaf->at(y_leaf->size() - 1);
            w_t z_w = y_leaf->weight_at(y_leaf->size() - 1);
            y_leaf->remove(y_leaf->size() - 1);
            lf->insert(0, z_bit, z_w);
            child_i++;
          } else {
            bool z_bit = y_leaf->at(0);
            w_t z_w = y_leaf->weight_at(0);
            y_leaf->remove(0);
            lf->push_back(z_bit, z_w);
          }
        } else {
          if (y_is_prev) {
            for (size_t k = 0; k < y_leaf->size(); ++k) {
              lf->insert(k, y_leaf->at(k), y_leaf->weight_at(k));
            }
            child_i += y_leaf->size();
            delete y_leaf;
            x->leaves.erase(x->leaves.begin() + j - 1);
            x->nr_children_--;
            j--;
          } else {
            for (size_t k = 0; k < y_leaf->size(); ++k) {
              lf->push_back(y_leaf->at(k), y_leaf->weight_at(k));
            }
            delete y_leaf;
            x->leaves.erase(x->leaves.begin() + j + 1);
            x->nr_children_--;
          }
        }
      }

      x->leaves[j]->remove(child_i);
      x->rebuild(0);

      node* tmp_child = x;
      node* tmp_parent = x->parent_;

      while (tmp_parent != nullptr) {
        tmp_parent->rebuild(0);
        tmp_child = tmp_parent;
        tmp_parent = tmp_child->parent_;
      }
    } else {
      node* new_ch_root = x->children[j]->remove(child_i);

      if (new_ch_root != nullptr) {
        delete x->children[j];

        x->children[j] = new_ch_root;
        new_ch_root->set_parent(x);
        new_ch_root->set_rank(j);
      }

      x->rebuild(0);
    }

    node* new_root = nullptr;
    if (this->is_root() && !this->has_leaves_ && this->nr_children_ == 1) {
      new_root = this->children[0];
      new_root->set_parent(nullptr);

      this->children.clear();
      this->nr_children_ = 0;
    }

    return new_root;
  }

  void free_mem() {
    if (has_leaves_) {
      for (uint32_t i = 0; i < nr_children_; ++i) delete leaves[i];

    } else {
      for (uint32_t i = 0; i < nr_children_; ++i) {
        children[i]->free_mem();
        delete children[i];
      }
    }

    leaves.clear();
    children.clear();

    nr_children_ = 0;
  }

  void set_rank(uint32_t r) { rank_ = r; }
  void set_parent(node* new_p) { parent_ = new_p; }

  size_t size_in_bits() const {
    size_t bits = sizeof(*this) * 8;
    bits += children.capacity() * sizeof(node*) * 8;
    bits += leaves.capacity() * sizeof(leaf_type*) * 8;
    if (has_leaves_) {
      for (size_t i = 0; i < nr_children_; ++i) {
        if (leaves[i]) {
          bits += leaves[i]->size_in_bits();
        }
      }
    } else {
      for (size_t i = 0; i < nr_children_; ++i) {
        if (children[i]) {
          bits += children[i]->size_in_bits();
        }
      }
    }
    return bits;
  }

  size_t size_in_bytes() const {
    return (size_in_bits() + 7) / 8;
  }

private:
  static uint64_t free_capacity(const leaf_type& l) {
    if (l.size() >= 2 * B_leaf) return 0;
    return 2 * B_leaf - l.size();
  }

  bool leaf_can_loose(const leaf_type* l) const {
    return l->size() >= B_leaf + 1;
  }

  bool is_root() const { return parent_ == nullptr; }

  bool can_loose_child() const {
    return is_root() || (nr_children_ >= B_fan_out + 2);
  }

  bool is_full() const {
    return nr_children_ == MAX_CHILDREN;
  }

  /*
   * helper functions for child search
   */
  uint64_t find_child(uint64_t i) const {
    uint64_t j = 0;
    while (j + 1 < nr_children_ && subtree_sizes[j] <= i) {
      j++;
    }
    return j;
  }

  uint64_t find_1(uint64_t x) const {
    uint64_t j = 0;
    while (j + 1 < nr_children_ && subtree_psums[j] < x) {
      j++;
    }
    return j;
  }

  uint64_t find_0(uint64_t x) const {
    uint64_t j = 0;
    while (j + 1 < nr_children_ && (subtree_sizes[j] - subtree_psums[j]) < x) {
      j++;
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

  void new_children(uint32_t i, node* left, node* right) {
    assert(i < nr_children_);
    assert(!is_full());

    children.insert(children.begin() + i + 1, right);
    children[i] = left;

    nr_children_++;

    for (uint32_t k = i; k < nr_children_; ++k) {
      children[k]->set_rank(k);
      children[k]->set_parent(this);
    }

    rebuild(0);
  }

  void new_children(uint32_t i, leaf_type* left, leaf_type* right) {
    assert(i < nr_children_);
    assert(!is_full());

    leaves.insert(leaves.begin() + i + 1, right);
    leaves[i] = left;

    nr_children_++;

    rebuild(0);
  }

  void insert_without_split(uint64_t i, bool bit, w_t weight) {
    assert(i <= size());

    uint32_t j = (i == size() ? nr_children_ - 1 : find_child(i));
    uint64_t prev_sz = (j == 0 ? 0 : subtree_sizes[j - 1]);
    uint64_t pos = i - prev_sz;

    if (has_leaves_) {
      if (free_capacity(*leaves[j]) == 0) {
        leaf_type* right_leaf = leaves[j]->split();

        if (pos <= leaves[j]->size()) {
          leaves[j]->insert(pos, bit, weight);
        } else {
          right_leaf->insert(pos - leaves[j]->size(), bit, weight);
        }

        new_children(j, leaves[j], right_leaf);
      } else {
        leaves[j]->insert(pos, bit, weight);

        rebuild(j);
      }
    } else {
      node* new_ch_root = children[j]->insert(pos, bit, weight);

      if (new_ch_root != nullptr) {
        children[j] = new_ch_root;

        new_ch_root->set_parent(this);
        new_ch_root->set_rank(j);
      }

      rebuild(j);
    }
  }

  node* split() {
    assert(is_full());

    uint32_t mid = nr_children_ / 2;
    node* right = nullptr;

    if (has_leaves_) {
      std::vector<leaf_type*> r_leaves(leaves.begin() + mid, leaves.end());
      leaves.erase(leaves.begin() + mid, leaves.end());

      right = new node(std::move(r_leaves), parent_, rank_ + 1);
    } else {
      std::vector<node*> r_children(children.begin() + mid, children.end());
      children.erase(children.begin() + mid, children.end());

      right = new node(std::move(r_children), parent_, rank_ + 1);
    }

    nr_children_ = mid;

    rebuild(0);

    return right;
  }

  std::array<uint64_t, MAX_CHILDREN> subtree_sizes{};
  std::array<uint64_t, MAX_CHILDREN> subtree_psums{};
  std::array<w_t, MAX_CHILDREN> subtree_aggs{};

  std::vector<node*> children;
  std::vector<leaf_type*> leaves; // empty if vertex isn't leaf

  node* parent_ = nullptr;

  uint32_t rank_ = 0;
  uint32_t nr_children_ = 0;
  bool has_leaves_ = false;
};

#endif // W_SPSI_NODE_H