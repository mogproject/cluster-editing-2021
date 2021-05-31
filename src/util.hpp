#pragma once
#include <assert.h>
#include <signal.h>
#include <math.h>
#include <cstring>
#include <iostream>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <queue>
#include <stack>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <chrono>
#include <functional>
#include <random>

// Constants
int const INF = 1000000000;  // 1e9
int const B = 64;            // number of bits in one integer
int const MAX_N = 700;       // max number of vertices

// Packing integers
constexpr inline int to_pair(int a, int b) { return a * MAX_N + b; }
constexpr inline int from_pair_0(int p) { return p / MAX_N; }
constexpr inline int from_pair_1(int p) { return p % MAX_N; }

constexpr inline int to_triple(int a, int b, int c) {
  return (a < 0 || b < 0 || c < 0) ? -1 : (a * MAX_N + b) * MAX_N + c;
}
constexpr inline int from_triple_0(int t) { return t < 0 ? -1 : t / MAX_N / MAX_N; }
constexpr inline int from_triple_1(int t) { return t < 0 ? -1 : t / MAX_N % MAX_N; }
constexpr inline int from_triple_2(int t) { return t < 0 ? -1 : t % MAX_N; }
