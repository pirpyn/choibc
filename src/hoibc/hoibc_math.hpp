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

  template<typename T>
  matrix<T> operator*(const matrix<T>& A, const T& x){
    matrix<T> C;
    C[0][0] = A[0][0]*x;
    C[1][0] = A[1][0]*x;
    C[0][1] = A[0][1]*x;
    C[1][1] = A[1][1]*x;
    return C;
  }

  template<typename T>
  matrix<T> operator*(const T& x, const matrix<T>& A){
    return A*x;
  }

  template<typename T>
  matrix<T> operator/(const matrix<T>& A, const T& x){
    matrix<T> C;
    C[0][0] = A[0][0]/x;
    C[1][0] = A[1][0]/x;
    C[0][1] = A[0][1]/x;
    C[1][1] = A[1][1]/x;
    return C;
  }

  template<typename T>
  matrix<T> operator+(const matrix<T>& A, const matrix<T>& B){
    matrix<T> C;
    C[0][0] = A[0][0] + B[0][0];
    C[1][0] = A[1][0] + B[1][0];
    C[0][1] = A[0][1] + B[0][1];
    C[1][1] = A[1][1] + B[1][1];
    return C;
  }

  template<typename T>
  matrix<T> operator-(const matrix<T>& A, const matrix<T>& B){
    matrix<T> C;
    C[0][0] = A[0][0] - B[0][0];
    C[1][0] = A[1][0] - B[1][0];
    C[0][1] = A[0][1] - B[0][1];
    C[1][1] = A[1][1] - B[1][1];
    return C;
  }

  template<typename T>
  matrix<T> matmul(const matrix<T>& A, const matrix<T>& B){
    matrix<T> C;
    C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
    C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
    C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
    C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
    return C;
  }

  template<typename T>
  matrix<T> inv(const matrix<T>& A){
    matrix<T> B;
    T delta = A[0][0]*A[1][1] - A[1][0]*A[0][1];

    B[1][1] = 1./delta * A[0][0];
    B[0][1] = -1./delta * A[0][1];
    B[1][0] = -1./delta * A[1][0];
    B[0][0] = 1./delta * A[1][1];

    return B;
  }

}
#endif