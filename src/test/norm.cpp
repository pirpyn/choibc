#include <iostream>
#include "../hoibc/hoibc_math.hpp"
#include <cmath>
#define CHECK(expr) if (expr) { cout << "SUCCESS" << endl;} else { cerr << "FAIL" << endl; no_of_errors++;};

using namespace std;


int main(){
  cout << "Testing norm" << endl;
  int no_of_errors {0};

  // TODO floating point number should not be compared exactly.
  CHECK(hoibc::norm(hoibc::array<hoibc::real> {})==0.)
  CHECK(hoibc::norm(hoibc::array<hoibc::integer> {3,4})==5)
  CHECK(hoibc::norm(hoibc::array<hoibc::real> {0.,0.,1.})==1.)
  CHECK(hoibc::norm(hoibc::array<hoibc::real> {1.,1.,1.})==std::sqrt(3.))
  CHECK(hoibc::norm(hoibc::array<hoibc::complex> {hoibc::complex(0.,0.),hoibc::complex(0.,1.)})==1.)
  CHECK(hoibc::norm(hoibc::matrix<hoibc::real> {{ {{1.,0.}} , {{0.,1.}} }})==std::sqrt(2.))
  CHECK(hoibc::norm(hoibc::matrix<hoibc::complex> {{ {{hoibc::complex(1.,0.),hoibc::complex(0.,0.)}} , {{hoibc::complex(0.,0.),hoibc::complex(0.,1.)}} }})==std::sqrt(2.))
  CHECK(hoibc::norm(hoibc::matrix<hoibc::complex> {{ {{hoibc::complex(1.,0.),hoibc::complex(1.,0.)}} , {{hoibc::complex(0.,0.),hoibc::complex(0.,1.)}} }})==std::sqrt(3.))

  return no_of_errors;
}