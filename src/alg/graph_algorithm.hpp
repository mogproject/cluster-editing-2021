#pragma once

#include "../data/component_info.hpp"

namespace mog {
namespace alg {

class GraphAlgorithm {
 public:
  template <int L>
  static std::vector<mog::data::ComponentInfo<L>> get_component_info(mog::data::Graph<L> const& g,
                                                                     typename mog::data::Graph<L>::BM mask) {
    std::vector<mog::data::ComponentInfo<L>> ret;
    while (!mask.empty()) {
      auto root = mask.front();
      auto ci = bfs(g, root);
      ret.push_back(ci);
      mask -= ci.mask;
    }
    return std::move(ret);
  }

  template <int L>
  static mog::data::ComponentInfo<L> bfs(mog::data::Graph<L> const& g, int root) {
    mog::data::ComponentInfo<L> ci;

    std::queue<int> q;
    q.push(root);
    ci.mask |= root;
    ci.n = 1;

    while (!q.empty()) {
      auto v = q.front();
      q.pop();
      // auto odd = oddity[v];

      BITMAP_FOREACH_START(g.nbr1[v], u)
      ++ci.m;
      if (ci.mask[u]) {
        // already visited
      } else {
        // new node
        ++ci.n;
        ci.mask |= u;
        q.push(u);
      }
      BITMAP_FOREACH_END;
    }

    ci.m /= 2;
    return std::move(ci);
  }

  /**
   * Complexity: O(n^3)?
   */
  template <int L>
  static int pack_p3(mog::data::Graph<L> const& g) {
    int ret = 0;

    typename mog::data::Graph<L>::SM avail;

    // (0) Setup: O(n^2 L)
    auto nbr2 = g.closed_neighborhood_r2() - g.nbr1 - mog::data::identity_matrix<L * B>;

    BITMAP_FOREACH_START(g.vertex_set, v)
    avail[v] |= (g.nbr1[v] | nbr2[v]) & g.vertex_set;
    BITMAP_FOREACH_END;
    assert(avail.is_symmetric());

    // (1) Process locked edges: O(n^2 L)
    BITMAP_FOREACH_START(g.vertex_set, v)
    auto us = (g.nbr1[v] - g.perm[v]) & avail[v];  // v's locked neighbors
    BITMAP_FOREACH_START(us, u)
    auto ws1 = ((g.nbr1[v] - g.nbr1[u] - u) & avail[v]) & avail[u];
    ret += ws1.size();

    BITMAP_FOREACH_START(ws1, w)
    if (!(g.perm[v][w] || g.perm[u][w])) return INF;  // infeasible
    if (g.perm[v][w]) avail -= std::make_pair(v, w);
    if (g.perm[u][w]) avail -= std::make_pair(u, w);
    BITMAP_FOREACH_END;

    auto ws2 = ((g.nbr1[u] - g.nbr1[v] - v) & avail[v]) & avail[u];
    ret += ws2.size();

    BITMAP_FOREACH_START(ws2, w)
    if (!(g.perm[v][w] || g.perm[u][w])) return INF;  // infeasible
    if (g.perm[v][w]) avail -= std::make_pair(v, w);
    if (g.perm[u][w]) avail -= std::make_pair(u, w);
    BITMAP_FOREACH_END;

    BITMAP_FOREACH_END;
    BITMAP_FOREACH_END;

    // (2) Process locked non-edges
    BITMAP_FOREACH_START(g.vertex_set, v)
    auto us = (avail[v] - g.nbr1[v]) - g.perm[v];  // v's locked non-neighbors
    BITMAP_FOREACH_START(us, u)
    auto ws = ((g.nbr1[v] & g.nbr1[u]) & avail[v]) & avail[u];  // v and u's common neighbors
    ret += ws.size();

    BITMAP_FOREACH_START(ws, w)
    if (!(g.perm[v][w] || g.perm[u][w])) return INF;  // infeasible
    if (g.perm[v][w]) avail -= std::make_pair(v, w);
    if (g.perm[u][w]) avail -= std::make_pair(u, w);

    BITMAP_FOREACH_END;

    BITMAP_FOREACH_END;
    BITMAP_FOREACH_END;

    // (3) Process others
    int num_avail = 0;
    std::vector<int> active_degrees(mog::data::Graph<L>::N);
    typename mog::data::Graph<L>::BM candidates_set;

    // compute active vertices
    BITMAP_FOREACH_START(g.vertex_set, v)
    if (!avail[v].empty()) {
      num_avail += avail[v].size();
      candidates_set |= v;
    }
    BITMAP_FOREACH_END;

    num_avail /= 2;

    // compute active degrees
    BITMAP_FOREACH_START(candidates_set, v)
    auto degv = (g.nbr1[v] & candidates_set).size();  // count only active neighbors

    auto p = to_pair(degv, v);
    active_degrees[v] = p;
    BITMAP_FOREACH_END;

    while (num_avail >= 3 && !candidates_set.empty()) {
      int p = INF;
      BITMAP_FOREACH_START(candidates_set, vv)
      p = std::min(p, active_degrees[vv]);
      BITMAP_FOREACH_END;
      assert(p != INF);
      int v = from_pair_1(p);  // choose a vertex with the min-degree

      auto ws = g.nbr1[v] & avail[v];
      BITMAP_FOREACH_START(ws, w)
      auto u = ((g.nbr1[w] & avail[w] & avail[v]) - g.nbr1[v]).front();
      if (u >= 0) {
        // found P_3
        ++ret;
        avail ^= std::make_pair(v, u);
        avail ^= std::make_pair(u, w);
        avail ^= std::make_pair(v, w);
        num_avail -= 3;
      }
      BITMAP_FOREACH_END;

      assert(avail.is_symmetric());

      // no valid P_3 including v => invalidate all edges indicent to v
      // Note: faster than `avail -= mog::data::Graph<L>::SM::cross_matrix(v);`
      num_avail -= avail[v].size();
      BITMAP_FOREACH_START(avail[v], u)
      avail[u] -= v;
      BITMAP_FOREACH_END;
      avail[v].clear();

      candidates_set -= v;

      // refresh information of v's active neighbors
      BITMAP_FOREACH_START(g.nbr1[v] & candidates_set, uu)
      active_degrees[uu] -= to_pair(1, 0);       // decrement degree
      if (active_degrees[uu] < to_pair(1, 0)) {  // degree zero
        candidates_set -= uu;
      }
      BITMAP_FOREACH_END;
    }
    return ret;
  }

  template <int L>
  static std::vector<int> order_active_vertices(mog::data::Graph<L> const& g, typename mog::data::Graph<L>::BM const& candidates_set) {
    std::vector<int> ret;

    BITMAP_FOREACH_START(candidates_set, v)
    auto degv = (g.nbr1[v] & candidates_set).size();  // count only active neighbors
    ret.push_back(to_pair(degv, v));
    BITMAP_FOREACH_END;

    std::sort(ret.begin(), ret.end());
    return std::move(ret);
  }
};

}  // namespace alg
}  // namespace mog