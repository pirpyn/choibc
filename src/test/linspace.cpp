// #include <iostream>
// #include "../hoibc/hoibc_math.hpp"

// #define MACRO { cout << "SUCCESS" << endl;} else { cerr << "FAIL" << endl; no_of_errors++;};

// using namespace std;


// int main(){
//   cout << "Testing linspace" << endl;
//   int no_of_errors {0};
//   if ((hoibc::array<hoibc::real> {})==hoibc::linspace(0.,1.,0)) MACRO
//   if ((hoibc::array<hoibc::real> {0.})==hoibc::linspace(0.,1.,1)) MACRO
//   if ((hoibc::array<hoibc::real> {0.,1.})==hoibc::linspace(0.,1.,2)) MACRO
//   if ((hoibc::array<hoibc::real> {0.,.5,1.})==hoibc::linspace(0.,1.,3)) MACRO
//   if ((hoibc::array<hoibc::real> {1.,3.})==hoibc::linspace(1,3,2)) MACRO
//   if ((hoibc::array<hoibc::real> {0.,-1.})==hoibc::linspace(0,-1.,2)) MACRO

//   return no_of_errors;

// }