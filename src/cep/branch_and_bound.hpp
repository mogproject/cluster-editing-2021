#pragma once

#include "local_search.hpp"
#include "search_context.hpp"
#include "../data/progress.hpp"

namespace mog {
namespace cep {

/**
 * Branch and Bound.
 */
class BranchAndBound {
 public:
  template <int L>
  static std::pair<bool, SearchResult<L>> run(SearchContext<L> &ctx, int known_best) {
    typedef mog::data::Graph<L> G;

    int best = known_best;
    typename G::SM cert;  // labels are shuffled

    long long node_limit = ctx.node_count + ctx.node_limit;
    int max_queue_size = ctx.search_queue.size();

    while (!ctx.search_queue.empty()) {
      //------------------------------------------------------
      // Phase 1: Metaprocessing
      //------------------------------------------------------

      max_queue_size = std::max(max_queue_size, (int)ctx.search_queue.size());

      ++ctx.node_count;
      if (ctx.node_count >= node_limit) { break; }

      //------------------------------------------------------
      // Phase 2: Choosing Next Node
      //------------------------------------------------------

      auto &item = ctx.search_queue.top();
      auto &g = item.g;
      ctx.search_queue.pop();

      // Immediate return
      if (g.num_edited >= best) {
        ctx.progress.add(item.weight);
        continue;
      }

      //------------------------------------------------------
      // Phase 3: Reduction
      //------------------------------------------------------

      auto red_ret = Reducer::reduce_and_lock(g, false, false);
      if (!red_ret) {
        // infeasible instance
        ctx.progress.add(item.weight);
        continue;
      }

      //------------------------------------------------------
      // Phase 4: Found a solution
      //------------------------------------------------------
      if (g.vertex_set.empty()) {
        assert(G({g.n, g.nbr1.to_vector()}).is_cluster_graph());  // sanity check

        if (g.num_edited < best) {
          // update best
          best = g.num_edited;
          cert = g.nbr1;
        }
        ctx.progress.add(item.weight);
        continue;
      }

      //------------------------------------------------------
      // Phase 5: Compute a lowerbound
      //------------------------------------------------------
      int lb = 0;
      if (best != INF) { lb = mog::alg::GraphAlgorithm::pack_p3(g); }
      if (g.num_edited + lb >= best) {
        ctx.progress.add(item.weight);
        continue;
      }

      //------------------------------------------------------
      // Phase 6: Branch
      //------------------------------------------------------
      ctx.brancher.run2(g, best, item.depth, item.last_branch_on, item.last_branch_index, item.weight, ctx.search_queue);
    }

    //------------------------------------------------------
    // Convert result
    //------------------------------------------------------
    bool finished = ctx.search_queue.empty();
    SearchResult<L> ret;
    ret.first = best;

    if (best < known_best) { ret.second = cert.map_indices(ctx.reverse_map); }
    return std::make_pair(finished, ret);
  }

 private:
};
}  // namespace cep
}  // namespace mog