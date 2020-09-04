#include <iostream>
#include "../hoibc/hoibc_math.hpp"
#include <cmath>
#include <limits>

#define CHECK(label) if (abs(computed-reference) < 10.*std::numeric_limits<double>::epsilon()) { \
        cout << "[OK]   " << label << ", error = " << abs(computed - reference) << endl; \
    } else { \
        cerr << "[FAIL] " << label << ", error = " << abs(computed - reference) << endl; \
        no_of_errors++;\
    };

using namespace std;
using namespace hoibc;

int main(){
  cout << "Testing norm" << endl;
  int no_of_errors {0};

  const hoibc::complex c0 = hoibc::complex(0.,0.);
  const hoibc::complex c1 = hoibc::complex(1.,0.);
  const hoibc::complex ci = hoibc::complex(0.,1.);

  hoibc::real reference;
  hoibc::real computed;

  reference = 0.;
  computed  = norm(hoibc::array<hoibc::real> {});
  CHECK("[]           ");

  computed = norm(hoibc::array<hoibc::integer> {3,4});
  reference = 5;
  CHECK("[3,4]        ");

  computed  = norm(hoibc::array<hoibc::real> {0.,0.,1.});
  reference = 1.;
  CHECK("[0,0,1]      ");

  computed  = norm(hoibc::array<hoibc::real> {1.,1.,1.});
  reference = sqrt(3.);
  CHECK("[1,1,1]      ");

  computed  = norm(hoibc::array<hoibc::complex> {c0,ci});
  reference = 1.;
  CHECK("[0,i]        ");

  computed  = norm(matrix<hoibc::real> {{ {{1.,0.}} , {{0.,1.}} }});
  reference = sqrt(2.);
  CHECK("[1,i]        ");

  computed  = norm(matrix<hoibc::complex> {{ {{c1,c0}} , {{c0,ci}} }});
  reference = sqrt(2.);
  CHECK("[[1,0],[0,i]]");

  computed  = norm(matrix<hoibc::complex> {{ {{c1,c1}} , {{c0,ci}} }});
  reference = sqrt(3.);
  CHECK("[[1,1],[0,i]]");

  return no_of_errors;
}