#pragma once

#include "graph.hpp"

namespace mog {
namespace data {

template <int L>
class ComponentInfo {
  typedef mog::data::Graph<L> G;

 public:
  typename G::BM mask;
  // std::vector<mog::data::Edge> bridges;
  int n;
  int m;
  // bool is_bipartite;

  // ComponentInfo() : n(0), m(0), is_bipartite(true) {}
  ComponentInfo() : n(0), m(0) {}

  bool is_clique() const { return m * 2 == n * (n - 1); }
};
}  // namespace data
}  // namespace mog