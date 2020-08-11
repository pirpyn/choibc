#include "hoibc_ibc0.hpp"

using namespace hoibc;

void hoibc::hoibc_ibc0::get_coeff_no_suc(const std::vector<real>& f1, const std::vector<real>& f2, const big_matrix<complex>& gex, const real& k0){
  switch (this->mode)  {
  case 2:
    this->coeff["a0"] = gex[0][0][0][0];
    break;
  default:
    std::cerr << "hoibc_ibc0::get_coeff_no_suc: mode = " << this->mode << " unknown" << std::endl;
    break;
  }
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

void hoibc::hoibc_ibc0::disp_coeff(std::ostream& out){
  out << "# Z = a0*I" << std::endl;
  print_complex(out,this->coeff["a0"],"    a0");
}