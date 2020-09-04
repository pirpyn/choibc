#include <iostream>
#include "../hoibc/hoibc_bessel.hpp"
#include "../hoibc/hoibc_types.hpp"
#include <limits>

#define CHECK(label) if (abs(computed-reference) < 10.*std::numeric_limits<double>::epsilon()) { \
        cout << "[OK]   " << label << ", error = " << abs(computed - reference) << endl; \
    } else { \
        cerr << "[FAIL] " << label << ", error = " << abs(computed - reference) << endl; \
        no_of_errors++;\
    }

using namespace std;
using namespace hoibc;

int main(){
  cout << "Testing bessel1" << endl;
  int no_of_errors {0};
  hoibc::complex reference;
  hoibc::complex computed;
  
  reference = 0.765197686557966551449717526102663220909274289755325241861;
  computed  = bessel1(0.,hoibc::complex(1.,0.));
  CHECK("bessel1(0.,1.)   ");

  reference = 0.019985850304223122424228390950848990680633578859027929558;
  computed  = bessel1(0.,hoibc::complex(100.,0.));
  CHECK("bessel1(0.,100.) ");

  reference = 2.6306151236874532069978536877905029440885704143207273E-10;
  computed  = bessel1(10.,hoibc::complex(1.,0.));
  CHECK("bessel1(10.,1.)  ");

  reference = 0.096366673295861559674314024870401848311755419825021855917;
  computed  = bessel1(100.,hoibc::complex(100.,0.));
  CHECK("bessel1(10.,100.)");

  reference = -0.440050585744933515959682203718914913;
  computed  = bessel1p(0.,hoibc::complex(1.,0.));
  CHECK("bessel1p(0.,1.)");

  return no_of_errors;

}