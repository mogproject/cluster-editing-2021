#pragma once
#include "../io/parser.hpp"
#include "matrix.hpp"

#define ASSERT_RANGE_2D(x, y) assert(0 <= (x) && (x) < N && 0 <= (y) && (y) < N)

namespace mog {
namespace data {

typedef std::pair<int, int> Edge;

/**
 * Graph representation.
 */
template <int L>
class Graph {
 public:
  static constexpr int const N = L * B;
  typedef Matrix<N, true> SM;  // symmetric matrix
  typedef Bitmap<N> BM;

  int n;                 // number of vertices, including both active and inactive ones
  int num_edited;        // number of edits so far
  Bitmap<N> vertex_set;  // active vertices

  Matrix<N, true> nbr1;  // adjacency matrix (symmetric); may include inactive vertices
  // Matrix<N, true> nbr2;  // radius-2 neighbors (symmetric); may include inactive vertices
  Matrix<N, true> perm;  // editable vertex pairs; may NOT include inactive vertices

 public:
  Graph() {}

  /**
   * Construct a graph from input.
   *
   * Complexity: O(n^2 L)
   */
  Graph(mog::io::EdgeListInput const& el) : num_edited(0), n(el.n) {
    // edges
    for (auto& e : el.edges) nbr1 |= e;

    // vertices (compute degrees)
    for (int v = 0; v < el.n; ++v) vertex_set |= v;

    // r2 neighborhood
    auto nbr2 = closed_neighborhood_r2() - nbr1 - identity_matrix<N>;

    // permission
    perm = nbr1 | nbr2;
  }

  /**
   * Returns the number of active vertices.
   *
   * Complexity: O(nL)
   */
  size_t num_vertices() const { return vertex_set.size(); }

  /**
   * Complexity: O(nL)
   */
  size_t num_edges() const { return nbr1.size(); }

  /**
   * Complexity: O(n^2 L)
   */
  size_t num_triangles() const {
    size_t ret = 0;
    for (int i = 0; i < n - 1; ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (nbr1[i][j]) { ret += (nbr1[i] & nbr1[j]).size(); }
      }
    }
    return ret / 3;
  }

  size_t num_active_triangles() const {
    size_t ret = 0;
    BITMAP_FOREACH_START(vertex_set, i)
    BITMAP_FOREACH_START(vertex_set, j)
    if (i < j) {
      if (nbr1[i][j]) { ret += (nbr1[i] & nbr1[j]).size(); }
    }
    BITMAP_FOREACH_END;
    BITMAP_FOREACH_END;
    return ret / 3;
  }

  /**
   * Complexity: O(n^2 L)
   */
  double triangle_ratio() const {
    if (n <= 2) return 0.0;
    return num_triangles() * 6.0 / n / (n - 1) / (n - 2);
  }

  double active_triangle_ratio() const {
    int nn = num_vertices();
    if (nn <= 2) return 0.0;
    return num_active_triangles() * 6.0 / nn / (nn - 1) / (nn - 2);
  }

  /**
   * Complexity: O(L)
   */
  size_t num_common_neighbors(int v, int u) const { return (nbr1[v] & nbr1[u]).size(); }

  /**
   * Complexity: O(L)
   */
  size_t num_noncommon_neighbors(int v, int u) const { return ((nbr1[v] ^ nbr1[u]) - v - u).size(); }

  /**
   * This includes inactive vertices.
   *
   * Complexity: O(L)
   */
  int degree(int v) const { return nbr1[v].size(); }

  int active_degree(int v) const { return (nbr1[v] & vertex_set).size(); }

  /**
   * Complexity: O(L)
   */
  std::vector<int> get_vertices() const {
    std::vector<int> ret;
    BITMAP_FOREACH_START(vertex_set, v)
    ret.push_back(v);
    BITMAP_FOREACH_END;
    return std::move(ret);
  }
  /**
   * Complexity: O(L)
   */
  BM get_vertices_mask() const { return vertex_set; }

  /**
   * Complexity: O(n^2)
   */
  std::vector<mog::data::Edge> get_edges() const { return nbr1.to_vector(); }
  /**
   * Complexity:
   */
  std::vector<BM> components() const { return components(vertex_set, nbr1); }

  static std::vector<BM> components(BM const& vs, SM const& nbr) {
    std::vector<BM> ret;
    BM visited;

    BITMAP_FOREACH_START(vs, v)
    if (visited[v]) continue;

    // BFS
    BM c = BM(v), q(v);
    visited |= v;

    while (true) {
      auto u = q.pop_front();
      if (u == -1) break;

      auto frontier = nbr[u] - visited;  // may contain vertices not in vertex_set
      if (frontier.empty()) continue;

      q |= frontier;
      c |= frontier;
      visited |= frontier;
    }
    ret.push_back(c);
    BITMAP_FOREACH_END;
    return std::move(ret);
  }

  bool is_clique(BM const& vs) const {
    int cnt = 0;

    BITMAP_FOREACH_START(vs, v)
    cnt += degree(v);
    BITMAP_FOREACH_END;
    return cnt == vs.size() * (vs.size() - 1);
  }

  bool is_cluster_graph() const {
    auto cc = components();
    return std::all_of(cc.begin(), cc.end(), [this](BM const& c) { return this->is_clique(c); });
  }

  Graph<L> relabel_vertices(std::vector<int> const& mapping) const {
    assert(mapping.size() == n);
    std::vector<mog::data::Edge> edges;
    for (auto& p : nbr1.to_vector()) { edges.push_back(std::make_pair(mapping[p.first], mapping[p.second])); }
    return Graph<L>(mog::io::EdgeListInput({n, edges}));
  }

  //------------------------------------------------------
  // Graph Modification
  //------------------------------------------------------

 private:
  /**
   * Complexity: O(?)
   */
  bool handle_lock(int v_, int u_, bool has_edge_) {
    assert(v_ != u_);

    std::queue<int> q;
    q.push(to_triple(has_edge_ ? 1 : 0, v_, u_));

    while (!q.empty()) {
      int x = q.front();
      q.pop();

      bool has_edge = from_triple_0(x) != 0;
      int v = from_triple_1(x);
      int u = from_triple_2(x);

      // (1) Already locked
      if (!perm[v][u]) {
        if (has_edge == nbr1[v][u]) continue;  // do nothing
        return false;                          // infeasible
      }

      // (2) Edge modification
      if (has_edge != nbr1[v][u]) {
        ++num_edited;
        nbr1 ^= std::make_pair(v, u);

        if (has_edge) {
          // (2-a) Edge addition
          // nbr2 ^= std::make_pair(v, u);

          // if u and w <- N[v] are not adjacent, then w-u becomes distance 2
          // for (int i = 0; i < 2; ++i) {
          //   auto a = i == 0 ? v : u;
          //   auto b = i == 0 ? u : v;
          //   // BITMAP_FOREACH_START(nbr1[a] - nbr1[b] - nbr2[b] - b, w)
          //   // assert(!perm[w][b]);  // must have been locked as d(w,b) used to be >2
          //   // nbr2 |= std::make_pair(w, b);
          //   // BITMAP_FOREACH_END;
          // }
        } else {
          // (2-b) Edge deletion
          // if (!(nbr1[v] & nbr1[u]).empty()) nbr2 |= std::make_pair(v, u);  // d_G(v,u)=1 -> d_G'(v,u)=2

          for (int i = 0; i < 2; ++i) {
            auto a = i == 0 ? v : u;
            auto b = i == 0 ? u : v;

            // if w-v-u is the only shortest path, w-u becomes a permanent non-edge
            // BITMAP_FOREACH_START(nbr1[a] & nbr2[b], w)
            BITMAP_FOREACH_START(nbr1[a] - nbr1[b], w)
            if ((nbr1[w] & nbr1[b]).empty()) {
              // nbr2 -= std::make_pair(w, b);
              q.push(to_triple(0, w, b));  // permanent non-edge
            }
            BITMAP_FOREACH_END;
          }
        }
      }

      // (3) Edge lock
      if (has_edge) {
        // (3-a) Edge lock
        assert(nbr1[v][u]);
        // fprintf(stderr, "edge lock: v=%d, u=%d\n", v, u);

        for (int i = 0; i < 2; ++i) {
          auto a = i == 0 ? v : u;
          auto b = i == 0 ? u : v;

          // w <- v's permanent neighbor & not u's permanent neighbor
          BITMAP_FOREACH_START((nbr1[a] - perm[a]) - (nbr1[b] - perm[b]) - b, w)
          if (!perm[b][w]) return false;  // permanent non-edge then infeasible
          // fprintf(stderr, "permanent edge bw=> v=%d, u=%d, a=%d, b=%d, w=%d\n", v, u, a, b, w);
          q.push(to_triple(1, b, w));  // otherwise, b-w is a permanent edge
          BITMAP_FOREACH_END;

          /*
           * w <- v's permanent non-neighbor & not u's permanent non-neighbor
           *   [~N(a) - P(a)] & ~[~N(b) - P(b)]
           * = [~N(a) & ~P(a)] & ~[~N(b) & ~P(b)]
           * = ~N(a) & ~P(a) & [N(b) | P(b)]
           * = [N(b) | P(b)] - N(a) - P(a)
           */
          BITMAP_FOREACH_START((nbr1[b] | perm[b]) - nbr1[a] - perm[a] - a, w)
          if (!perm[b][w]) return false;  // permanent edge then infeasible
          // fprintf(stderr, "permanent non-edge bw=> v=%d, u=%d, a=%d, b=%d, w=%d\n", v, u, a, b, w);
          q.push(to_triple(0, b, w));  // otherwise, b-w is a permanent non-edge
          BITMAP_FOREACH_END;
        }
      } else {
        // (3-b) Non-edge lock
        assert(!nbr1[v][u]);

        for (int i = 0; i < 2; ++i) {
          auto a = i == 0 ? v : u;
          auto b = i == 0 ? u : v;

          /*
           * w <- v's permanent neighbor & not u's permanent non-neighbor
           *   [N(a) - P(a)] & ~[~N(b) - P(b)]
           * = [N(a) & ~P(a)] & ~[~N(b) & ~P(b)]
           * = [N(a) & ~P(a)] & [N(b) | P(b)]
           * = [N(b) | P(b)] & N(a) - P(a)
           */
          BITMAP_FOREACH_START(((nbr1[b] | perm[b]) & nbr1[a]) - perm[a], w)
          if (!perm[b][w]) return false;  // permanent edge then infeasible
          q.push(to_triple(0, b, w));     // otherwise, b-w is a permanent non-edge
          BITMAP_FOREACH_END;
        }
      }
      perm -= std::make_pair(v, u);  // lock v-u
    }
    return true;
  }

 public:
  /**
   * Complexity: O(1)
   */
  bool lock(int v, int u) {
    assert(v != u);
    return handle_lock(v, u, nbr1[v][u]);
  }

  /**
   * Complexity: O(?)
   */
  bool edit_edge(int v, int u) {
    assert(v != u);
    return handle_lock(v, u, !nbr1[v][u]);
  }

  bool make_permanent_edge(int v, int u) {
    assert(v != u);
    return handle_lock(v, u, true);
  }

  bool make_permanent_non_edge(int v, int u) {
    assert(v != u);
    return handle_lock(v, u, false);
  }

  /**
   * Inactivates a vertex.
   *
   * Complexity: O(L)
   * Precondition: perm[v] must be empty; that is, all the edges incident to vertex v must be locked
   */
  inline void hide_vertex(int v) {
    assert(perm[v].empty());
    vertex_set -= v;  // deactivate vertex v
  }

  /**
   * @return true if a clique is made, false if the graph cannot make the given clique.
   */
  bool make_clique(BM const& vs) {
    auto cm = Matrix<N, true>::cross_matrix(vs);
    auto cl = Matrix<N, true>::clique_matrix(vs);
    int k = vs.size();  // clique size

    if (!((nbr1 & cm) - cl).empty()) return false;  // connected to outside

    auto after = nbr1 | cl;
    auto upd = nbr1 ^ after;
    if (!upd.subset(perm)) return false;  // no permission

    nbr1 |= cl;
    num_edited += upd.size();
    vertex_set -= vs;  // inactivate
    perm -= cm;        // non-editable
    return true;
  }

  std::string repr() const {
    std::stringstream ss;
    ss << "Graph(";
    ss << "n=" << n;
    ss << ", num_edited=" << num_edited;
    ss << ", nbr1=" << nbr1.repr();
    ss << ", perm=" << perm.repr();
    ss << ")";
    return ss.str();
  }

  /**
   * Complexity: O(n^2 L)
   */
  Matrix<N, true> closed_neighborhood_r2() const {
    Matrix<N, true> ret;
    BITMAP_FOREACH_START(vertex_set, v)
    BITMAP_FOREACH_START(nbr1[v], u)
    ret[v] |= nbr1[u];
    BITMAP_FOREACH_END;
    BITMAP_FOREACH_END;
    assert(ret.is_symmetric());
    return std::move(ret);
  }
};

}  // namespace data
}  // namespace mog