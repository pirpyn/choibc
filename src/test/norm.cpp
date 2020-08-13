#include <iostream>
#include "../hoibc/hoibc_math.hpp"

#define CHECK(expr) if (expr) { cout << "SUCCESS" << endl;} else { cerr << "FAIL" << endl; no_of_errors++;};

using namespace std;


int main(){
  cout << "Testing norm" << endl;
  int no_of_errors {0};
  CHECK(hoibc::norm(vector<hoibc::real> {})==0.)
  CHECK(hoibc::norm(vector<hoibc::integer> {3,4})==5)
  CHECK(hoibc::norm(vector<hoibc::real> {0.,0.,1.})==1.)
  CHECK(hoibc::norm(vector<hoibc::complex> {hoibc::complex(0.,0.),hoibc::complex(0.,1.)})==1.)

  return no_of_errors;

}