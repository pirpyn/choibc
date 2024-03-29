#include "hoibc_math.hpp"

namespace hoibc {

  matrix<complex> operator*(const matrix<real>& A, const complex& x){
    matrix<complex> C;
    C[0][0] = A[0][0]*x;
    C[1][0] = A[1][0]*x;
    C[0][1] = A[0][1]*x;
    C[1][1] = A[1][1]*x;
    return C;
  }

  matrix<complex> operator*(const complex& x, const matrix<real>& A){
    return A*x;
  }

  big_matrix<complex> operator*(const big_matrix<real>& A, const complex& x){
    std::size_t n1 = A.size();
    std::size_t n2 = A[0].size();
    big_matrix<complex> C = big_init(n1,n2,complex(0.,0.));
    for (std::size_t i=0;i<n1;i++){
      for (std::size_t j=0;j<n2;j++){
        C[i][j] = A[i][j]*x;
      }
    }
    return C;
  }

  big_matrix<complex> operator*(const complex& x, const big_matrix<real>& A){
    return A*x;
  }

  matrix<complex> conj(const matrix<complex>& A){
    matrix<complex> B = A;
    B[0][0] = std::conj(A[0][0]);
    B[0][1] = std::conj(A[0][1]);
    B[1][0] = std::conj(A[1][0]);
    B[1][1] = std::conj(A[1][1]);
    return B;
  }

  matrix<complex> operator+(const matrix<complex>& A, const matrix<real>& B){
    matrix<complex> C;
    C[0][0] = A[0][0] + B[0][0];
    C[1][0] = A[1][0] + B[1][0];
    C[0][1] = A[0][1] + B[0][1];
    C[1][1] = A[1][1] + B[1][1];
    return C;
  }

  matrix<complex> operator+(const matrix<real>& A, const matrix<complex>& B){
    return B + A;
  }

  matrix<complex> operator*(const matrix<complex>& A, const matrix<real>& B){
    matrix<complex> C;
    C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
    C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
    C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
    C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
    return C;
  }

  matrix<complex> operator*(const matrix<real>& A, const matrix<complex>& B){
    return B * A;
  }

  big_matrix<complex> operator+(const big_matrix<complex>& A, const big_matrix<real>& B){
    big_matrix<complex> C = A;
    for (std::size_t i=0;i<A.size();i++){
      for (std::size_t j=0;j<A[i].size();j++){
        C[i][j] = A[i][j] + B[i][j];
      }
    }
    return C;
  }

  big_matrix<complex> operator+(const big_matrix<real>& A, const big_matrix<complex>& B){
    return B + A;
  }

  big_matrix<complex> operator*(const big_matrix<real>& A, const big_matrix<complex>& B){
    big_matrix<complex> C = B;
    for (std::size_t i=0;i<A.size();i++){
      for (std::size_t j=0;j<A[i].size();j++){
        C[i][j] = A[i][j] * B[i][j];
      }
    }
    return C;
  }

  big_matrix<complex> operator*(const big_matrix<complex>& A, const big_matrix<real>& B){
    return B*A;
  }
}