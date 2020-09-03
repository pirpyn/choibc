#include <iostream>
#include "../hoibc/hoibc_math.hpp"
#include "../hoibc/hoibc_types.hpp"

#define MACRO { cout << "SUCCESS" << endl;} else { cerr << "FAIL" << endl; no_of_errors++;}

using namespace std;
using namespace hoibc;

int main(){
  cout << "Testing overloading of * when multiplying a array<real> with scalar" << endl;
  int no_of_errors {0};
  const hoibc::array<hoibc::real> u = {1.,1.};
  if ((hoibc::array<hoibc::real> {1.,1.})==(u*1)) MACRO;
  if ((hoibc::array<hoibc::real> {2.,2.})==(u*2)) MACRO; 
  if ((hoibc::array<hoibc::real> {2.5,2.5})==(2.5*u)) MACRO; 
  if ((hoibc::array<hoibc::real> {-.5,-.5})==(-1*u*.5)) MACRO;

  return no_of_errors;

}