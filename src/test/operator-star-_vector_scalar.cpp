#include <iostream>
#include "../hoibc/hoibc_math.hpp"
#include "../hoibc/hoibc_types.hpp"

#define CHECK(label) if ((computed == reference).min()) { \
        cout << "[OK]   " << label << endl; \
    } else { \
        cerr << "[FAIL] " << label << endl; \
        no_of_errors++;\
    }
  
using namespace std;
using namespace hoibc;
int main(){
  cout << "Testing overloading of * when multiplying a array<real> with scalar" << endl;
  int no_of_errors {0};
  hoibc::array<hoibc::real> base, reference, computed;
  base = {1.,1.};

  computed = base*1;
  reference = {1.,1.};
  CHECK("{1.,1.}     == 1*{1.,1.}");

  computed = base*2;
  reference = {2.,2.};
  CHECK("{2.,2.}     == {1.,1.}*2");

  computed = 2.5*base;
  reference = {2.5,2.5};
  CHECK("{2.5.,2.5.} == 2.5*{1.,1.}");

  computed = -1*base*.5;
  reference = {-.5,-.5};
  CHECK("{-.5.,-.5}  == -1*{1.,1.}*.5");
  return no_of_errors;

}