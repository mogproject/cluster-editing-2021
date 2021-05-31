#pragma once

#include "../data/graph.hpp"

namespace mog {
namespace cep {

template <int L>
class Brancher {
  typedef mog::data::Graph<L> G;

  G const *root_;
  std::default_random_engine gen_;
  int strategy_;
  int num_samples_;
  bool sampling_order_;                                  // true if high-degree vertices should be chosen
  std::uniform_real_distribution<double> epsilon_dist_;  // used for random tie-breaking of metrics
  std::uniform_real_distribution<double> pref_dist_;     // used for random tie-breaking of preference
  bool cert_provided_;
  typename G::SM cert_;
  typename G::SM removed_edges_;
  typename G::SM added_edges_;
  std::vector<unsigned int> ordered_edges_;
  std::map<unsigned int, int> ordered_edges_rev_;

 public:
  /**
   * Constructor.
   */
  Brancher() {}

  Brancher(uint32_t seed,              //
           int strategy,               //
           int num_samples,            //
           G const *root,              //
           typename G::SM const *cert  //
           )
      : gen_(seed),
        root_(root),
        strategy_(strategy),
        num_samples_(num_samples),
        epsilon_dist_(0, 1e-6),
        pref_dist_(-0.5, 0.5) {
    sampling_order_ = strategy < 10;

    // compute edited edges
    cert_provided_ = false;
    if (cert) {
      cert_provided_ = true;
      cert_ = *cert;
      removed_edges_ = root->nbr1 - cert_;
      added_edges_ = cert_ - root->nbr1;
    }

    // precompute ordering
    if (strategy >= 10) {
      assert(cert_provided_);
      order_all_edges();
    }
  }

  void order_all_edges() {
    if (strategy_ % 10 >= 5) {
      // component-wise
      auto cc = G::components(root_->vertex_set, cert_);

      // sort by component size (largest first)
      std::sort(cc.begin(), cc.end(),
                [](typename G::BM const &a, typename G::BM const &b) { return a.size() > b.size(); });
      int cindex = 0;
      for (auto &c : cc) {
        BITMAP_FOREACH_START(c, i)
        for (int j = i + 1, jj = 0; j < root_->n; ++j, ++jj) {
          if (!root_->perm[i][j]) continue;
          int is_external = c[j] ? 0 : 1;
          // ordered_edges_.push_back(to_triple(to_pair(is_external, cindex), i, j));  // internal first, component-wise
          // ordered_edges_.push_back(to_triple(to_pair(is_external, jj), i, j));  // internal first, spanning
          ordered_edges_.push_back(to_triple(to_pair(1 - is_external, cindex), i, j));  // external first, component-wise
        }
        BITMAP_FOREACH_END;
        ++cindex;
      }
    } else if (strategy_ % 10 < 5) {
      // ordering
      unsigned int order_non_edge_edited = 0;
      unsigned int order_non_edge_non_edited = 1;
      unsigned int order_edge_edited = 0;
      unsigned int order_edge_non_edited = 1;

      if (strategy_ % 10 == 1) {
        order_non_edge_edited = 1;
        order_non_edge_non_edited = 0;
        order_edge_edited = 1;
        order_edge_non_edited = 0;
      }
      bool high_degree_first = true;

      for (int i = 0; i < root_->n; ++i) {
        auto degi = root_->degree(i);
        // auto degi = cert[i].size();
        for (int j = i + 1; j < root_->n; ++j) {
          if (!root_->perm[i][j]) continue;

          auto degsum = root_->degree(j) + degi;
          // auto degsum = degi;
          int weight = high_degree_first ? 2 * MAX_N - degsum : degsum;
          // int weight = high_degree_first ? 2 * MAX_N - degsum : degsum;

          if (removed_edges_[i][j]) {
            weight = to_pair(order_non_edge_edited * 2, weight);
          } else if (added_edges_[i][j]) {
            weight = to_pair(order_edge_edited * 2, weight);
          } else {
            auto bonus = (root_->nbr1[i][j] ? order_edge_non_edited : order_non_edge_non_edited) * 2;
            weight = to_pair(bonus, weight);
          }
          ordered_edges_.push_back(to_triple(weight, i, j));
        }
      }
    }

    std::sort(ordered_edges_.begin(), ordered_edges_.end());  // increasing order
    for (int i = 0; i < ordered_edges_.size(); ++i) {
      // reverse map
      auto val = ordered_edges_[i] % (MAX_N * MAX_N);
      ordered_edges_rev_[val] = i;
    }
  }

  bool prefer_non_edge(G const &g, int v, int u) {
    switch (strategy_ / 10) {
      case 1:
      case 5: return !cert_[v][u];  // solution state first
      case 2:
      case 6: return !g.nbr1[v][u];  // keep original state
      case 3:
      case 7: return cert_[v][u];  // solution state last
      case 4:
      case 8: return g.nbr1[v][u];  // edit first
      default: {
        if (g.num_common_neighbors(v, u) < g.num_noncommon_neighbors(v, u)) return true;
        if (g.num_common_neighbors(v, u) > g.num_noncommon_neighbors(v, u)) return false;
      }
    }
    return pref_dist_(gen_) < 0;
  }

  /**
   * Performs binary branching.
   */
  template <typename T>
  void run(G const &g, int known_best, std::stack<T> &search_queue) {
    // choose target edge
    auto target = choose_target(g);
    if (target < 0) return;

    // create next states
    int v = from_pair_0(target), u = from_pair_1(target);
    bool can_have_edge = g.num_edited + (g.nbr1[v][u] ? 0 : 1) + g.num_noncommon_neighbors(v, u) < known_best;
    bool can_have_nonedge = g.num_edited + (g.nbr1[v][u] ? 1 : 0) + g.num_common_neighbors(v, u) < known_best;

    if (can_have_edge && can_have_nonedge) {
      G h1(g), h2(g);
      h1.make_permanent_edge(v, u);
      h2.make_permanent_non_edge(v, u);
      if (prefer_non_edge(g, v, u)) {
        search_queue.push(h1);  // edge
        search_queue.push(h2);  // nonedge
      } else {
        search_queue.push(h2);  // nonedge
        search_queue.push(h1);  // edge
      }
    } else if (can_have_edge) {
      G h(g);
      h.make_permanent_edge(v, u);
      search_queue.push(h);
    } else if (can_have_nonedge) {
      G h(g);
      h.make_permanent_non_edge(v, u);
      search_queue.push(h);
    }
  }

  template <typename T>
  void run2(G const &g, int known_best, int depth, int last_branch_on, int last_branch_index, int weight, std::stack<T> &search_queue) {
    // choose target edge
    bool preordered = strategy_ >= 50;
    bool use_preordered = preordered;
    int target = -1;
    if (preordered) {
      // preordered
      target = choose_preordered_edge(g, last_branch_index + 1);
    } else if (strategy_ >= 10) {
      // hybrid
      int last_v = from_triple_1(last_branch_on);
      int last_u = from_triple_2(last_branch_on);
      if (last_branch_on >= 0 && root_->nbr1[last_v][last_u] == g.nbr1[last_v][last_u]) {
        // choose an edge that causes at least one edge modification
        target = choose_prefered_edge(g, last_v, last_u);
      }
      if (target < 0) {
        target = choose_preordered_edge(g, last_branch_index + 1);
        use_preordered = true;
      }
    } else {
      // dynamic
      target = choose_target(g, last_branch_on);
    }

    if (target < 0) return;

    // create next states
    int v = use_preordered ? from_triple_1(ordered_edges_[target]) : from_pair_0(target);
    int u = use_preordered ? from_triple_2(ordered_edges_[target]) : from_pair_1(target);

    bool can_have_edge = g.num_edited + (g.nbr1[v][u] ? 0 : 1) + g.num_noncommon_neighbors(v, u) < known_best;
    bool can_have_nonedge = g.num_edited + (g.nbr1[v][u] ? 1 : 0) + g.num_common_neighbors(v, u) < known_best;
    auto nxt_bi = use_preordered ? target : last_branch_index;

    if (can_have_edge && can_have_nonedge) {
      G h1(g), h2(g);
      h1.make_permanent_edge(v, u);
      h2.make_permanent_non_edge(v, u);
      if (prefer_non_edge(g, v, u)) {
        search_queue.push({depth + 1, to_pair(v, u), nxt_bi, weight + 1, h1});  // edge
        search_queue.push({depth + 1, to_pair(v, u), nxt_bi, weight + 1, h2});  // nonedge
      } else {
        search_queue.push({depth + 1, to_pair(v, u), nxt_bi, weight + 1, h2});  // nonedge
        search_queue.push({depth + 1, to_pair(v, u), nxt_bi, weight + 1, h1});  // edge
      }
    } else if (can_have_edge) {
      G h(g);
      h.make_permanent_edge(v, u);
      search_queue.push({depth + 1, to_pair(v, u), nxt_bi, weight, h});
    } else if (can_have_nonedge) {
      G h(g);
      h.make_permanent_non_edge(v, u);
      search_queue.push({depth + 1, to_pair(v, u), nxt_bi, weight, h});
    }
  }

 private:
  int choose_target(G const &g, int last_branch_on) {
    if (last_branch_on < 0) return choose_target(g);  // first branching

    int v = from_triple_1(last_branch_on);
    int u = from_triple_2(last_branch_on);
    bool has_edge = g.nbr1[v][u];

    if (root_->nbr1[v][u] != has_edge) return choose_target(g);  // edge was edited

    // limit the candidates
    std::vector<int> candidates;
    if (has_edge) {
      auto ws = g.nbr1[v] - g.perm[v] - g.nbr1[u];  // v's exclusive editable neighbor
      BITMAP_FOREACH_START(ws, w)
      candidates.push_back(to_pair(v, w));
      assert(v != w);
      BITMAP_FOREACH_END;

      auto ys = g.nbr1[u] - g.perm[u] - g.nbr1[v];  // u's exclusive editable neighbor
      BITMAP_FOREACH_START(ys, y)
      candidates.push_back(to_pair(u, y));
      assert(u != y);
      BITMAP_FOREACH_END;
    } else {
      auto ws = (g.nbr1[v] - g.perm[v]) & (g.nbr1[u] - g.perm[u]);  // editable common neighbor
      BITMAP_FOREACH_START(ws, w)
      candidates.push_back(to_pair(v, w));
      candidates.push_back(to_pair(u, w));
      BITMAP_FOREACH_END;
    }

    if (candidates.empty()) return choose_target(g);  // no such candidates

    // find the best edge
    int best_edge = -1;
    double best_metric = -INF;

    for (auto x : candidates) {
      auto metric = compute_metric(g, from_pair_0(x), from_pair_1(x));
      if (metric > best_metric) {
        best_edge = x;
        best_metric = metric;
      }
    }
    assert(best_edge >= 0);
    return best_edge;
  }

  double compute_metric(G const &g, int v, int u) {
    double metric;
    switch (strategy_) {
      case 0: {
        int common = g.num_common_neighbors(v, u);
        int noncom = g.num_noncommon_neighbors(v, u);
        int noncomx = noncom + 1 - (g.nbr1[v][u] ? 1 : 0);
        metric = std::abs(common - noncomx);
        break;
      }
      // case 1: {
      //   // Strategy 1: measure the degree sum
      //   metric = from_pair_0(candidate_vertices[i]) + g.active_degree(u) + epsilon_dist_(gen_);
      //   break;
      // }
      case 3: {
        // Strategy 3: measure the number of common neighbors
        metric = g.num_common_neighbors(v, u) + epsilon_dist_(gen_);
        break;
      }
      case 4: {
        // Strategy 4 measure the number of noncommon neighbors
        metric = g.num_noncommon_neighbors(v, u) + epsilon_dist_(gen_);
        break;
      }
      case 6: {
        // Strategy 6 edited edges first
        unsigned int const order_non_edge_edited = 0;
        unsigned int const order_non_edge_non_edited = 1;
        unsigned int const order_edge_edited = 0;
        unsigned int const order_edge_non_edited = 1;

        int bonus = order_non_edge_non_edited;  // non-edge non-edited
        if (removed_edges_[v][u]) {
          bonus = order_non_edge_edited;  // non-edge edited
        } else if (added_edges_[v][u]) {
          bonus = order_edge_edited;  // edge edited
        } else if (g.nbr1[v][u]) {
          bonus = order_edge_non_edited;  // edge non-edited
        }
        metric = (4 - bonus) * 1000 + g.num_common_neighbors(v, u) + epsilon_dist_(gen_);
        break;
      }
      case 7: {
        // Strategy 7 edges first
        unsigned int const order_non_edge_edited = 0;
        unsigned int const order_non_edge_non_edited = 0;
        unsigned int const order_edge_edited = 1;
        unsigned int const order_edge_non_edited = 1;

        int bonus = order_non_edge_non_edited;  // non-edge non-edited
        if (removed_edges_[v][u]) {
          bonus = order_non_edge_edited;  // non-edge edited
        } else if (added_edges_[v][u]) {
          bonus = order_edge_edited;  // edge edited
        } else if (g.nbr1[v][u]) {
          bonus = order_edge_non_edited;  // edge non-edited
        }
        metric = (4 - bonus) * 1000 + epsilon_dist_(gen_);
        break;
      }
      default: assert(false);
    }
    return metric;
  }

  int choose_target(G const &g) {
    // collect active degree of all vertices
    std::vector<int> candidate_vertices;

    BITMAP_FOREACH_START(g.vertex_set, v)
    if (!g.perm[v].empty()) { candidate_vertices.push_back(to_pair(g.active_degree(v), v)); }
    BITMAP_FOREACH_END;

    // sort vertices by degree
    if (sampling_order_) {
      // largest first
      std::sort(candidate_vertices.begin(), candidate_vertices.end(), std::greater<int>());
      assert(candidate_vertices.size() < 2 || from_pair_0(candidate_vertices[0]) >= from_pair_0(candidate_vertices[1]));
    } else {
      std::sort(candidate_vertices.begin(), candidate_vertices.end(), std::less<int>());
    }

    // find the best edge
    int best_edge = -1;
    double best_metric = -INF;

    int num_samples = std::min(num_samples_, (int)candidate_vertices.size());
    for (int i = 0; i < num_samples; ++i) {
      int v = from_pair_1(candidate_vertices[i]);

      BITMAP_FOREACH_START(g.perm[v], u)
      auto metric = compute_metric(g, v, u);
      if (metric > best_metric) {
        best_edge = to_pair(v, u);
        best_metric = metric;
      }
      BITMAP_FOREACH_END;
    }

    return best_edge;
  }

  int choose_preordered_edge(G const &g, int index) {
    for (; index < ordered_edges_.size(); ++index) {
      auto x = ordered_edges_[index];
      if (!g.perm[from_triple_1(x)][from_triple_2(x)]) continue;  // not editable
      break;
    }
    return index >= ordered_edges_.size() ? -1 : index;
  }

  int choose_prefered_edge(G const &g, int last_v, int last_u) {
    std::vector<int> candidates;

    auto f = [&](int a, int b) { candidates.push_back(to_pair(std::min(a, b), std::max(a, b))); };

    if (g.nbr1[last_v][last_u]) {
      // edge between v and u => w-v-u or v-u-w is a P_3
      auto ws = g.nbr1[last_v] - g.perm[last_v] - g.nbr1[last_u];
      BITMAP_FOREACH_START(ws, w)
      f(last_v, w);
      BITMAP_FOREACH_END;
      auto ys = g.nbr1[last_u] - g.perm[last_u] - g.nbr1[last_v];
      BITMAP_FOREACH_START(ys, w)
      f(last_u, w);
      BITMAP_FOREACH_END;
    } else {
      // no edge between v and u => v-w-u is a P_3
      auto ws = (g.nbr1[last_v] - g.perm[last_v]) & (g.nbr1[last_u] - g.perm[last_u]);
      BITMAP_FOREACH_START(ws, w)
      f(last_v, w);
      f(last_u, w);
      BITMAP_FOREACH_END;
    }
    if (candidates.empty()) return -1;

    // find the earliest one in the preordered list
    int ret = -1;
    int best = INF;

    for (auto &e : candidates) {
      assert(ordered_edges_rev_.find(e) != ordered_edges_rev_.end());
      if (ordered_edges_rev_[e] < best) {
        best = ordered_edges_rev_[e];
        ret = e;
      }
    }

    return ret;
  }
};
}  // namespace cep
}  // namespace mog