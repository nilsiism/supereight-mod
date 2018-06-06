#ifndef SE_MATH_HELPER_H
#define SE_MATH_HELPER_H

#include "se_common.h"

namespace se {
  namespace math {
    template <typename T>
      T fracf(const T& v) {
        return v - v.array().floor().matrix();
      }

    template <typename T>
      T floorf(const T& v) {
        return v.array().floor();
      }

    template <typename T>
    inline T fabs(const T& v) {
        return v.cwiseAbs();
      }

    template <typename Scalar>
      inline Scalar sq(Scalar a) {
        return a*a;
      };

    template <typename Scalar>
    inline bool in(const Scalar v, const Scalar a, 
          const Scalar b) {
        return v >= a && v <= b;
      }

    constexpr int log2_const(int n){
      return (n < 2 ? 0 : 1 + log2_const(n/2));
    }
  }
}

#endif
