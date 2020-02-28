#ifndef _H_HOIBC_MATH
#define _H_HOIBC_MATH

#include <vector>
#include "hoibc_types.hpp"

namespace hoibc {

  template<typename Tstart, typename Tend> 
  std::vector<real> linspace(const Tstart& start, const Tend& end, const int& n) {

    std::vector<real> v;
    if (n<1) return v;
    v.resize(n);
    v[0] = start;
    if (n<2) {
      return v;
    }
    else {
      for ( unsigned int i = 0; i < v.size(); i++) {
        v[i] = static_cast<real>(start) + i*(static_cast<real>(end)-static_cast<real>(start))/(n-1);
      }
    }
    return v;
  }

  template<typename T>
  std::vector<real> operator*(const std::vector<real>& v, const T& x){
    std::vector<real> w (v);
    for (auto& y : w) {
      y *= x;
    }
    return w;
  }

  template<typename T>
  std::vector<real> operator*(const T& x, const std::vector<real>& v){
      return v*x;
  }

  template<typename T>
  std::vector<real> operator/(const std::vector<real>& v, const T& x){
    std::vector<real> w (v);
    for (auto& y : w) {
      y /= x;
    }
    return w;
  }

}
#endif