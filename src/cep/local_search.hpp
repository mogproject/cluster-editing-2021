#pragma once

#include "../data/graph.hpp"
#include "cep.hpp"

namespace mog {
namespace cep {

template <int L>
class LocalSearcher {
 private:
  typedef mog::data::Graph<L> G;
  G const *root_;

 public:
  LocalSearcher(G const *root) : root_(root) {}

  SearchResult<L> run(typename G::SM edges) {
    bool success = true;
    int best_sofar = (root_->nbr1 ^ edges).size();

    while (success) {
      success = false;

      for (int i = 0; i < root_->n; ++i) {
        auto nbri = edges[i];
        auto nbri_sz = edges[i].size();  // size of i's clique minus one
        auto cnbri = edges[i] | i;       // all of i's clique members

        //------------------------------------------------------
        // 1-Separation
        //------------------------------------------------------
        if (nbri_sz > 1) {  // clique size must be at least 3
          auto sx = (root_->nbr1[i] & nbri).size();
          if (2 * sx <= nbri_sz) {
            // prefer isolation to inclusion
            success = true;
            // isolate vertex i
            edges -= G::SM::cross_matrix(i);
            // update best if it is strictly better
            if (2 * sx < nbri_sz) { best_sofar = (root_->nbr1 ^ edges).size(); }
          }
        }

        //------------------------------------------------------
        // 1-Transfer
        //------------------------------------------------------
        typename G::BM chosen = cnbri;

        for (int j = 0; j < root_->n; ++j) {
          if (chosen[j]) continue;

          auto cnbrj = edges[j] | j;
          chosen |= cnbrj;

          auto trynbr = edges;
          trynbr -= G::SM::cross_matrix(i);

          BITMAP_FOREACH_START(cnbrj, x)
          trynbr |= std::make_pair(x, i);
          BITMAP_FOREACH_END;

          auto sz = (root_->nbr1 ^ trynbr).size();
          if (sz < best_sofar) {
            success = true;
            best_sofar = sz;
            edges = trynbr;
          }
        }
      }
    }

    return std::make_pair(best_sofar, edges);
  }
};
}  // namespace cep
}  // namespace mog