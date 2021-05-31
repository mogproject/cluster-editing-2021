#pragma once

#include "../data/graph.hpp"

namespace mog {
namespace cep {

/**
 * Search node.
 */
template <int L>
struct SearchNode {
  int depth = 0;
  int last_branch_on = -1;
  int last_branch_index = -1;
  int weight = 0;  // If this node is a leaf, then add 2^{-weight} to the overall progress.
  mog::data::Graph<L> g;
};

/**
 * Search result.
 */
template <int L>
using SearchResult = std::pair<int, typename mog::data::Graph<L>::SM>;

}  // namespace cep
}  // namespace mog