#pragma once

#include "../util.hpp"

namespace mog {
namespace data {
/**
 * Represents a precise real value between 0 and 1, inclusive,
 * generated by a sum of 2^{-k}, where 0 <= k < 250000.
 */
class Progress {
  static int const PRECISION = 4000;
  unsigned long long data_[PRECISION];

 public:
  Progress() { memset(data_, 0, sizeof(data_)); }

  double value() const { return data_[0] / pow(2.0, 63); }
  /**
   * Add 2^{-k} to the value.
   */
  void add(int k) {
    int y = k / 64;
    assert(y < PRECISION);

    unsigned long long carry = 1ull << (63 - k % 64);
    while (y >= 0 && carry > 0) {
      auto old = data_[y];
      data_[y] += carry;
      if (old <= data_[y]) break;
      carry = 1;
      --y;
    }
  }
};
}  // namespace data
}  // namespace mog