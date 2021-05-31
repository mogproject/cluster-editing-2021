#pragma once

#include "cep.hpp"
#include "branch.hpp"
#include "../data/progress.hpp"

namespace mog {
namespace cep {

/**
 * Search result.
 */
template <int L>
using SearchResult = std::pair<int, typename mog::data::Graph<L>::SM>;

/**
 * Search context.
 */
template <int L>
struct SearchContext {
  typedef mog::data::Graph<L> G;

  int id;
  G root;
  std::vector<int> label_map, reverse_map;
  int strategy;
  long long node_count;
  long long node_limit;
  mog::data::Progress progress;
  Brancher<L> brancher;
  std::stack<SearchNode<L>> search_queue;
  std::vector<unsigned int> ordered_edges;

  SearchContext(int id,                            //
                G const& root_,                    //
                uint32_t shuffle_seed,             //
                uint32_t branch_seed,              //
                int strategy,                      //
                int num_branch_samples,            //
                long long node_limit,              //
                SearchResult<L> const* known_best  //
                )
      : id(id), strategy(strategy), label_map(root_.n), reverse_map(root_.n), node_count(0), node_limit(node_limit) {
    // Shuffle vertices
    shuffle_vertices(root_, shuffle_seed);

    // Set up brancher
    if (known_best) {
      auto cert = known_best->second.map_indices(label_map);
      brancher = Brancher<L>(branch_seed, strategy, num_branch_samples, &root, &cert);
    } else {
      brancher = Brancher<L>(branch_seed, strategy, num_branch_samples, &root, nullptr);
    }

    // Push root node
    search_queue.push({0, -1, -1, 0, root});
  }

 private:
  void shuffle_vertices(G const& root_, uint32_t shuffle_seed) {
    for (int i = 0; i < root_.n; ++i) {
      label_map[i] = i;
      reverse_map[i] = i;
    }

    std::default_random_engine shuffle_gen(shuffle_seed);
    std::shuffle(label_map.begin(), label_map.end(), shuffle_gen);
    for (int i = 0; i < root_.n; ++i) reverse_map[label_map[i]] = i;  // update reverse map
    root = root_.relabel_vertices(label_map);
  }
};
}  // namespace cep
}  // namespace mog