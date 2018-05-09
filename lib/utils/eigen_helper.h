#ifndef EIGEN_HELPER_H
#define EIGEN_HELPER_H

#include <cmath>
#include "Eigen/Dense"

namespace octlib {
  namespace math {
    template <typename T>
      T fracf(const T& v) {
        return v - floor(v);
      }

    template <typename T>
      T floorf(const T& v) {
        return floor(v);
      }

    template <typename S, typename T>
      auto max(const S& a, const T& b) {
        return a.cwiseMax(b);
      }

    template <typename S, typename T>
      auto min(const S& a, const T& b) {
        return a.cwiseMin(b);
      }

    template <typename T>
      T fabs(const T& v) {
        return v.cwiseAbs();
      }

    template <typename T, typename U, typename V>
      T in(const T&v, const U& a, const V& b) {
        return (v.array() >= a.array()) && (v.array() <= b.array());
      }

    constexpr int log2_const(int n){
      return (n < 2 ? 0 : 1 + log2_const(n/2));
    }
  }
}

#endif
