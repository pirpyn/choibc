#include "hoibc_math.hpp"

using namespace hoibc;

hoibc::matrix<complex> hoibc::operator*(const matrix<real>& A, const complex& x){
  matrix<complex> C;
  C[0][0] = A[0][0]*x;
  C[1][0] = A[1][0]*x;
  C[0][1] = A[0][1]*x;
  C[1][1] = A[1][1]*x;
  return C;
}

hoibc::matrix<complex> hoibc::operator*(const complex& x, const matrix<real>& A){
  return A*x;
}

hoibc::big_matrix<complex> hoibc::operator*(const big_matrix<real>& A, const complex& x){
  std::size_t n1 = A.size();
  std::size_t n2 = A[0].size();
  big_matrix<complex> C = big_init(n1,n2,static_cast<complex>(0));
  for (std::size_t i=0;i<n1;i++){
    for (std::size_t j=0;j<n2;j++){
      C[i][j] = A[i][j]*x;
    }
  }
  return C;
}

hoibc::big_matrix<complex> hoibc::operator*(const complex& x, const big_matrix<real>& A){
  return A*x;
}
