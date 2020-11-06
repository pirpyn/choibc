#include <iostream>
#include "../hoibc/hoibc_math.hpp"

#define CHECK(label) \
    if ((computed.size()==reference.size() && computed.size() == 0) || (computed == reference).min()) { \
        cout << "[OK]   " << label << endl; \
    } else { \
        cerr << "[FAIL] " << label << endl; \
        no_of_errors++;\
    }

using namespace std;
using namespace hoibc;

int main(){
  cout << "Testing linspace" << endl;
  int no_of_errors {0};
  hoibc::array<hoibc::real> reference, computed;

  computed = linspace(0.,1.,0);
  reference = hoibc::array<hoibc::real> {};
  CHECK("{}         == linspace(0.,1.,0.)");

  computed = linspace(0.,1.,2);
  reference = hoibc::array<hoibc::real> {0.,1.};
  CHECK("{0.,1.}    == linspace(0.,1.,2)");

  computed = linspace(0.,1.,3);
  reference = hoibc::array<hoibc::real> {0.,.5,1.};
  CHECK("{0.,.5,1.} == linspace(0.,1.,3)");

  computed = linspace(1,3,2);
  reference = hoibc::array<hoibc::real> {1.,3.};
  CHECK("{1.,3.}    == linspace(1,3,2)");

  computed = linspace(0.,-1.,2);
  reference = hoibc::array<hoibc::real> {0.,-1.};
  CHECK("{0.,-1.}   == linspace(0.,-1.,2)");

  return no_of_errors;
}