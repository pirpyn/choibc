#include "hoibc_ibc0.hpp"

using namespace hoibc;

struct hoibc_ibc0::coeff {
  complex a0 = complex(0.,0.);
};

void hoibc::hoibc_ibc0::get_coeff_no_suc(const std::vector<real>& f1, const std::vector<real>& f2, const big_matrices<complex>& gex, const real& k0){
  std::cout << "hoibc_ibc0::get_coeff_no_suc: i do nothing" << std::endl;
}

void hoibc::hoibc_ibc0::get_impedance(){
  std::cout << "hoibc_ibc0::get_impedance: i do nothing" << std::endl;
}

void hoibc::hoibc_ibc0::array_to_coeff(){
  std::cout << "hoibc_ibc0::array_to_coeff: i do nothing" << std::endl;
}

void hoibc::hoibc_ibc0::coeff_to_array(){
  std::cout << "hoibc_ibc0::coeff_to_array: i do nothing" << std::endl;
}

void hoibc::hoibc_ibc0::get_suc(){
  std::cout << "hoibc_ibc0::get_suc: i do nothing" << std::endl;
}

void hoibc::hoibc_ibc0::disp_coeff(){
  std::cout << "hoibc_ibc0::disp_coeff: i do nothing" << std::endl;
}