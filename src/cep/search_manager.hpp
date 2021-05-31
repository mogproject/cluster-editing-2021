#pragma once

#include "rand_search.hpp"
#include "branch_and_bound.hpp"

namespace mog {
namespace cep {

template <int L>
class SearchManager {
 public:
  // Types
  typedef mog::data::Graph<L> G;
  static constexpr int const N = L * B;

  static std::vector<mog::data::Edge> run(G const& root, int seed = 12345) {
    // main logic
    //------------------------------------------------------
    // Round 1. Random Search
    //------------------------------------------------------
    std::default_random_engine gen(seed);
    std::uniform_int_distribution<uint32_t> int_dist({});

    SearchResult<L> best({INF, {}});

    bool large_graph = root.n >= 300;
    std::vector<int> strategies = {3, 4};
    int stability = 0;
    int const stability_threshold = large_graph ? 8 : 20;
    int allowance = INF;
    int num_branch_samples = 20;
    int num_max_solutions = 1;
    int num_tries = large_graph ? 10 : 120;
    uint32_t best_shuffle_seed = 0, best_branch_seed = 0;
    int best_strategy = strategies[0];
    bool finished = false;

    for (int i = 0; i < num_tries; ++i) {
      uint32_t shuffle_seed = int_dist(gen);
      uint32_t branch_seed = int_dist(gen);
      int strategy = strategies[i % strategies.size()];

      auto ret = RandSearch::run(root, shuffle_seed, branch_seed, num_max_solutions, strategy, num_branch_samples,
                                 allowance, nullptr, true);

      if (ret.second.first < best.first) {
        // update best
        best = ret.second;
        stability = 0;
        best_shuffle_seed = shuffle_seed;
        best_branch_seed = branch_seed;
        best_strategy = strategy;
      } else {
        if (++stability == stability_threshold) {
          if (num_max_solutions == 2) {
            // done with random search
            break;
          }
          num_max_solutions = 2;
          stability = 0;
        }
      }

      if (ret.first) {
        finished = true;
        break;  // complete
      }
    }

    //------------------------------------------------------
    // Round 2. Exhaustive Search
    //------------------------------------------------------

    if (!finished) {
      int num_contexts_per_strategy = 4;
      int initial_node_limit = 2000;
      int eliminate_thres = 50;

      if (root.n >= 300) {
        num_contexts_per_strategy = std::min(num_contexts_per_strategy, 1);
      } else if (root.n >= 200) {
        num_contexts_per_strategy = std::min(num_contexts_per_strategy, 2);
        eliminate_thres = 20;
      }

      // Create search contexts
      std::vector<SearchContext<L>> contexts;
      std::vector<double> progress_sofar;
      std::vector<int> bnb_strategies = {50, 51};

      for (auto st : bnb_strategies) {
        for (int j = 0; j < num_contexts_per_strategy; ++j) {
          auto ss = st == best_strategy && j == 0 ? best_shuffle_seed : int_dist(gen);
          auto bs = st == best_strategy && j == 0 ? best_branch_seed : int_dist(gen);
          auto ctx = SearchContext<L>({st * 100 + j, root, ss, bs, st, num_branch_samples, initial_node_limit, &best});
          contexts.push_back(ctx);
          progress_sofar.push_back(0);
        }
      }

      // main logic
      int index = 0;

      while (!finished) {
        auto bnb_result = BranchAndBound::run(contexts[index], best.first);
        if (bnb_result.second.first < best.first) {
          best = bnb_result.second;  // update best
        }
        if (bnb_result.first) {
          finished = true;  // finish
        }

        // context switch
        if (contexts.size() > 1) {
          ++index;
          if (index == contexts.size()) {
            if (true) {
              double diff = 0;

              // normalize limits
              for (int i = 0; i < contexts.size(); ++i) {
                auto& ctx = contexts[i];
                ctx.node_limit = std::min(100000LL, ctx.node_limit * 2);
                diff = std::max(diff, ctx.progress.value() - progress_sofar[i]);
              }

              // max progress diff bonus
              for (int i = 0; i < contexts.size(); ++i) {
                auto& ctx = contexts[i];
                double vv = ctx.progress.value();
                if (vv - progress_sofar[i] == diff) { ctx.node_limit *= 2; }
                progress_sofar[i] = vv;
              }

              // eliminate one context per strategy
              int jj = contexts.size() / strategies.size();
              if (jj > 1 && contexts[0].node_count > initial_node_limit * eliminate_thres) {
                for (int s = strategies.size() - 1; s >= 0; --s) {
                  double worst = 2;
                  int ww = -1;
                  for (int j = 0; j < jj; ++j) {
                    auto vv = contexts[s * jj + j].progress.value();
                    if (vv < worst) {
                      ww = s * jj + j;
                      worst = vv;
                    }
                  }
                  contexts.erase(contexts.begin() + ww);
                }
              }
            }
            index = 0;
          }
        }
      }
    }

    // collect result
    return (root.nbr1 ^ best.second).to_vector();
  }
};
}  // namespace cep
}  // namespace mog