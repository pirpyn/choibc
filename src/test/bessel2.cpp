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
  cout << "Testing bessel2" << endl;
  int no_of_errors {0};
  hoibc::complex reference;
  hoibc::complex computed;

  reference = hoibc::complex(0.7651976865579665514497175261026632209092742897553252418, -0.08825696421567695798292676602351516282781752309067554671);
  computed  = bessel2(0.,hoibc::complex(1.,0.));
  CHECK("bessel2(0.,1.)  ");

  reference = hoibc::complex(2.630615123687453206997853687791E-10,+1.2161801427868918928813042666797E8);
  computed  = bessel2(10.,hoibc::complex(1.,0.));
  CHECK("bessel2(10.,1.) ");

  reference = hoibc::complex(-0.4400505857449335159596822037189,-0.7812128213002887165471500000480);
  computed  = bessel2p(0.,hoibc::complex(1.,0.));
  CHECK("bessel2p(0.,1.) ");

  reference = hoibc::complex(2.618635056224421836032912981021E-9,-1.2093999378481599097712034732764E9);
  computed = bessel2p(10.,hoibc::complex(1.,0.));
  CHECK("bessel2p(10.,1.)");

  return no_of_errors;

}