#include "cep/search_manager.hpp"
#include "data/graph.hpp"
#include "io/parser.hpp"
#include "util.hpp"
#include <fstream>

using namespace std;
using namespace mog::io;
using namespace mog::data;
using namespace mog::cep;

EdgeListInput el;
std::vector<std::pair<int, int>> result;

template <int L>
void do_it() {
  int k = (el.n + B - 1) / B;
  if (k == L) {
    Graph<L> g(el);
    result = SearchManager<L>::run(g, RANDOM_SEED);
  }
}

int main(int argc, char *argv[]) {
  el = DIMACSParser::parse(cin);
  if (el.n == 0) return 0;  // do nothing

  if (el.n >= MAX_N) {
    fprintf(stderr, "too large n\n");
    return 1;
  }

  do_it<1>();
  do_it<2>();
  do_it<3>();
  do_it<4>();
  do_it<5>();
  do_it<6>();
  do_it<7>();
  do_it<8>();
  do_it<9>();
  do_it<10>();
  do_it<11>();

  for (auto &e : result) { printf("%d %d\n", e.first + 1, e.second + 1); }
  return 0;
}