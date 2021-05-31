#pragma once

#include "reduce.hpp"
#include "branch.hpp"
#include "local_search.hpp"

namespace mog {
namespace cep {

/**
 * Random search.
 */
class RandSearch {
 public:
  template <int L>
  static std::pair<bool, SearchResult<L>> run(mog::data::Graph<L> const& root, uint32_t shuffle_seed, uint32_t branch_seed,
                                              int num_max_solutions, int strategy, int num_branch_samples, int allowance,
                                              SearchResult<L> const* known_best, bool local_search_enabled) {
    typedef mog::data::Graph<L> G;

    //------------------------------------------------------
    // 1. Shuffle vertices
    //------------------------------------------------------
    std::vector<int> label_map(root.n), reverse_map(root.n);
    for (int i = 0; i < root.n; ++i) {
      label_map[i] = i;
      reverse_map[i] = i;
    }

    std::default_random_engine shuffle_gen(shuffle_seed);
    std::shuffle(label_map.begin(), label_map.end(), shuffle_gen);
    for (int i = 0; i < root.n; ++i) reverse_map[label_map[i]] = i;  // update reverse map
    auto root_shuffled = root.relabel_vertices(label_map);

    //------------------------------------------------------
    // 2. Set up Brancher
    //------------------------------------------------------
    Brancher<L> brancher(branch_seed, strategy, num_branch_samples, &root, nullptr);

    //------------------------------------------------------
    // 3. Branch and bound
    //------------------------------------------------------
    LocalSearcher<L> ls(&root_shuffled);

    int solution_count = 0;
    SearchResult<L> best({known_best ? known_best->first : INF, {}});  // reset the best value

    std::stack<G> search_queue;
    search_queue.push(root_shuffled);
    long long node_count = 0;

    while (!search_queue.empty() && solution_count < num_max_solutions) {
      ++node_count;

      auto g = search_queue.top();
      search_queue.pop();

      // Immediate return
      if (g.num_edited >= best.first + allowance) { continue; }

      // Reduce
      auto red_ret = Reducer::reduce_and_lock(g, true, true);
      if (!red_ret) {
        // infeasible instance
        continue;
      }

      if (!g.vertex_set.empty()) {
        // Computing lowerbound
        int lb = 0;
        if (best.first != INF) { lb = mog::alg::GraphAlgorithm::pack_p3(g); }

        // Branch
        if (g.num_edited + lb < best.first + allowance) {
          brancher.run(g, best.first + allowance, search_queue);
        } else {
          // branch cut
        }
      } else {
        // Found a solution
        ++solution_count;

        // sanity check
        assert(G({g.n, g.nbr1.to_vector()}).is_cluster_graph());

        // local search
        if (local_search_enabled) {
          auto ls_ret = ls.run(g.nbr1);
          // update best
          if (ls_ret.first < best.first) { best = ls_ret; }
        } else if (g.num_edited < best.first) {
          // update best
          best.first = g.num_edited;
          best.second = g.nbr1;
        }
      }
    }
    bool complete = search_queue.empty();

    //------------------------------------------------------
    // 4. Convert result
    //------------------------------------------------------
    SearchResult<L> ret = std::make_pair(best.first, best.second.map_indices(reverse_map));
    return std::make_pair(complete, ret);
  }

 private:
};
}  // namespace cep
}  // namespace mog