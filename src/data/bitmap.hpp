#pragma once
#include "../util.hpp"

namespace mog {
namespace data {

#define BITMAP_FOREACH_START(bm,j) {\
    for (int _i_ = 0; _i_ < L; ++_i_) { \
      if ((bm).data_[_i_]) { \
        auto _x_ = (bm).data_[_i_]; \
        while (_x_) { \
          auto (j) = _i_ * B + __builtin_ctzll(_x_); \
          _x_ &= _x_ - 1;

#define BITMAP_FOREACH_END }}}}

/**
 * Bitmap.
 */
template <int N>
class Bitmap {
 public:
  static constexpr int const L = (N + (B - 1)) / B;  // number of internal integers

  unsigned long long data_[L];
  

 private:
  inline void trim() {
    if (N % B) data_[L - 1] &= (1ULL << (N % B)) - 1;
  }

 public:
  Bitmap() { memset(data_, 0, sizeof(data_)); }
  /**
   * Constructs a singleton.
   */
  Bitmap(int x) : Bitmap() { *this |= x; }
  /**
   * Constructs from a list.
   */
  Bitmap(std::vector<int> const& xs) : Bitmap() {
    for (int x : xs) *this |= x;
  }

  inline void clear() {
    for (int i = 0; i < L; ++i) data_[i] = 0;
  }

  inline bool empty() const {
    for (int i = 0; i < L; ++i) {
      if (data_[i]) return false;
    }
    return true;
  }

  inline int front() const {
    int ret = -1;
    for (int i = 0; i < L; ++i) {
      if (data_[i]) {
        ret = i * B + __builtin_ctzll(data_[i]);
        break;
      }
    }
    return ret >= N ? -1 : ret;
  }

  inline int pop_front() {
    for (int i = 0; i < L; ++i)
      if (data_[i]) {
        int offset = __builtin_ctzll(data_[i]);
        int ret = i * B + offset;
        if (ret < N) {
          data_[i] ^= 1ULL << offset;
          return ret;
        }
      }
    return -1;
  }

  //--------------------------------------------------------
  //    Operators
  //--------------------------------------------------------

  /**
   * negation (not)
   */
  Bitmap<N> operator~() const {
    Bitmap<N> ret;
    for (int i = 0; i < L; ++i) ret.data_[i] = ~data_[i];
    ret.trim();
    return std::move(ret);
  }

  /**
   * set/union (or)
   */
  Bitmap<N>& operator|=(int x) {
    assert(0 <= x && x < N);
    data_[x / B] |= 1ULL << (x % B);
    return *this;
  }

  Bitmap<N>& operator|=(Bitmap<N> const& rhs) {
    for (int i = 0; i < L; ++i) data_[i] |= rhs.data_[i];
    return *this;
  }

  friend Bitmap<N> operator|(Bitmap<N> const& lhs, int x) {
    Bitmap<N> ret(lhs);
    ret |= x;
    return std::move(ret);
  }

  friend Bitmap<N> operator|(Bitmap<N> const& lhs, Bitmap<N> const& rhs) {
    Bitmap<N> ret;
    for (int i = 0; i < L; ++i) ret.data_[i] = lhs.data_[i] | rhs.data_[i];
    return std::move(ret);
  }

  /**
   * exclusive or (xor)
   */
  Bitmap<N>& operator^=(int x) {
    assert(0 <= x && x < N);
    data_[x / B] ^= 1ULL << (x % B);
    return *this;
  }

  Bitmap<N>& operator^=(Bitmap<N> const& rhs) {
    for (int i = 0; i < L; ++i) data_[i] ^= rhs.data_[i];
    return *this;
  }

  friend Bitmap<N> operator^(Bitmap<N> const& lhs, int x) {
    Bitmap<N> ret(lhs);
    ret ^= x;
    return std::move(ret);
  }

  friend Bitmap<N> operator^(Bitmap<N> const& lhs, Bitmap<N> const& rhs) {
    Bitmap<N> ret;
    for (int i = 0; i < L; ++i) ret.data_[i] = lhs.data_[i] ^ rhs.data_[i];
    return std::move(ret);
  }

  /**
   * intersection (and)
   */
  Bitmap<N> operator&=(int x) {
    assert(0 <= x && x < N);
    for (int i = 0; i < L; ++i) {
      if (i == x / B) {
        data_[x / B] &= 1ULL << (x % B);
      } else {
        data_[i] = 0;
      }
    }

    return *this;
  }

  Bitmap<N>& operator&=(Bitmap<N> const& rhs) {
    for (int i = 0; i < L; ++i) data_[i] &= rhs.data_[i];
    return *this;
  }

  friend Bitmap<N> operator&(Bitmap<N> const& lhs, int x) {
    Bitmap<N> ret(lhs);
    ret &= x;
    return std::move(ret);
  }

  friend Bitmap<N> operator&(Bitmap<N> const& lhs, Bitmap<N> const& rhs) {
    Bitmap<N> ret;
    for (int i = 0; i < L; ++i) ret.data_[i] = lhs.data_[i] & rhs.data_[i];
    return std::move(ret);
  }

  /**
   * reset/set minus
   */
  Bitmap<N>& operator-=(int x) {
    assert(0 <= x && x < N);
    data_[x / B] &= ~(1ULL << (x % B));
    return *this;
  }

  Bitmap<N>& operator-=(Bitmap<N> const& rhs) {
    for (int i = 0; i < L; ++i) data_[i] &= ~rhs.data_[i];
    return *this;
  }

  friend Bitmap<N> operator-(Bitmap<N> const& lhs, int x) {
    assert(0 <= x && x < N);
    Bitmap<N> ret(lhs);
    ret -= x;
    return std::move(ret);
  }

  friend Bitmap<N> operator-(Bitmap<N> const& lhs, Bitmap<N> const& rhs) {
    Bitmap<N> ret;
    for (int i = 0; i < L; ++i) ret.data_[i] = lhs.data_[i] & ~rhs.data_[i];
    return std::move(ret);
  }

  /**
   * get
   */
  bool operator[](int x) const {
    assert(0 <= x && x < N);
    return (data_[x / B] >> (x % B)) & 1ULL;
  }

  std::vector<int> to_vector() const {
    std::vector<int> ret;
    for (int i = 0; i < L; ++i) {
      if (data_[i]) {
        auto x = data_[i];
        while (x) {
          ret.push_back(i * B + __builtin_ctzll(x));
          x &= x - 1;
        }
      }
    }
    return std::move(ret);
  }

  /**
   * equality
   */
  friend inline bool operator==(Bitmap<N> const& lhs, Bitmap<N> const& rhs) {
    for (int i = 0; i < L; ++i) {
      if (lhs.data_[i] != rhs.data_[i]) return false;
    }
    return true;
  }

  friend inline bool operator!=(Bitmap<N> const& lhs, Bitmap<N> const& rhs) { return !(lhs == rhs); }

  int capacity() const { return N; }

  std::size_t size() const {
    std::size_t ret = 0;
    for (int i = 0; i < L; ++i) ret += __builtin_popcountll(data_[i]);
    return ret;
  }

  inline bool subset(Bitmap<N> const& rhs) const { return (*this & rhs) == *this; }

  inline bool superset(Bitmap<N> const& rhs) const { return rhs.subset(*this); }

  std::string encode() const {
    std::stringstream ss;
    for (int i = L - 1; i >= 0; --i) { ss << std::setfill('0') << std::setw(B / 4) << std::hex << data_[i]; }
    return ss.str();
  }
};

template <int N>
Bitmap<N> const empty_bitmap = Bitmap<N>();

template <int N>
std::ostream& operator<<(std::ostream& os, Bitmap<N> const& vm) {
  os << "Bitmap<" << N << ">({";
  int i = 0;
  for (auto x : vm.to_vector()) {
    if (i++ > 0) os << ",";
    os << x;
  }
  return os << "})";
}

}  // namespace data
}  // namespace mog