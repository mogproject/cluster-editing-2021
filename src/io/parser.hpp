#pragma once
#include <iostream>
#include <sstream>
#include <vector>

namespace mog {
namespace io {

class EdgeListInput {
 public:
  int n = 0;
  std::vector<std::pair<int, int>> edges;
};

class DIMACSParser {
 public:
  static EdgeListInput parse(std::istream &is) {
    EdgeListInput el;
    int m, v, u;

    for (std::string line; std::getline(is, line);) {
      if (line.empty()) continue;
      if (line[0] == 'c') continue;  // ignore comments

      auto ss = std::stringstream(line);

      if (line[0] == 'p') {
        std::string p, problem;

        ss >> p >> problem >> el.n >> m;
        continue;
      }

      // edge
      ss >> v >> u;
      el.edges.push_back({v - 1, u - 1});
    }
    return std::move(el);
  }
};

}  // namespace io
}  // namespace mog