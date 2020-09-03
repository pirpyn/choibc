#include <iostream>
#include "../hoibc/hoibc_bessel.hpp"
#include "../hoibc/hoibc_types.hpp"
#include <limits>

#define MACRO { cout << "SUCCESS" << endl;} else { cerr << "FAIL" << endl; no_of_errors++;}

using namespace std;
using namespace hoibc;

int main(){
  cout << "Testing bessel1 0 1" << endl;
  int no_of_errors {0};

  const hoibc::complex reference = 0.765197686557966551449717526102663220909274289755325241861;
  const hoibc::real error = abs(bessel1(0.,hoibc::complex(1.,0.))-reference);

  if (error > 10.*std::numeric_limits<double>::epsilon()) MACRO;

  return no_of_errors;

}