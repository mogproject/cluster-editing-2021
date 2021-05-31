#pragma once

#include "../alg/graph_algorithm.hpp"
#include "../data/component_info.hpp"
#include "cep.hpp"

namespace mog {
namespace cep {

class Reducer {
 public:
  /**
   * @return true if success, false if there is a contradiction
   *
   * Complexity: O(n^2 L)
   */
  template <int L>
  static bool reduce_and_lock(mog::data::Graph<L> &g, bool enable_twin_lock, bool enable_dense_reduction) {

    //------------------------------------------------------
    // Phase 1. True twin lock
    //------------------------------------------------------
    // An edge between true twins (u,v s.t. N[u]=N[v]) cannot be editable.
    //------------------------------------------------------
    if (enable_twin_lock) {
      BITMAP_FOREACH_START(g.vertex_set, v)
      auto cnbrv = g.nbr1[v] | v;
      BITMAP_FOREACH_START(g.nbr1[v] & g.perm[v], u)
      if (v < u && (cnbrv == (g.nbr1[u] | u))) {
        if (!g.perm[v][u]) continue;  // may have been updated during this loop
        // // sanity check
        assert(((g.nbr1[v] - g.perm[v]) & ((~g.nbr1[u]) - g.perm[u])).empty());
        assert(((g.nbr1[u] - g.perm[u]) & ((~g.nbr1[v]) - g.perm[v])).empty());

        auto lock_ret = g.lock(v, u);
        if (!lock_ret) {
          return false;  // infeasible
        }
      }
      BITMAP_FOREACH_END;
      BITMAP_FOREACH_END;
    }

    //------------------------------------------------------
    // Phase 3. Dense Component Reduction
    //------------------------------------------------------
    if (enable_dense_reduction) {
      std::queue<mog::data::ComponentInfo<L>> components;
      for (auto &c : mog::alg::GraphAlgorithm::get_component_info(g, g.vertex_set)) components.push(c);

      // process component by component
      while (!components.empty()) {
        auto c = components.front();
        components.pop();

        if (c.m * 2 >= (c.n - 1) * (c.n - 1)) {
          if (g.make_clique(c.mask)) {
            continue;  // done
          } else {
            // failure
          }
        }
      }
    }

    //------------------------------------------------------
    // Phase 4. Inactivate Vertices
    //------------------------------------------------------
    BITMAP_FOREACH_START(g.vertex_set, v)
    if (g.perm[v].empty()) { g.hide_vertex(v); }
    BITMAP_FOREACH_END;

    return true;
  }

 private:
};

}  // namespace cep
}  // namespace mog
