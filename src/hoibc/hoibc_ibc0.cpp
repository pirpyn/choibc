#include "hoibc_ibc0.hpp"
#include "hoibc_math.hpp"

using namespace hoibc;

void hoibc::hoibc_ibc0::get_coeff_no_suc(const std::vector<real>& f1, const std::vector<real>& f2, const big_matrix<complex>& gex, const real& k0){
  switch (this->mode)  {
  case 2:
    this->coeff.a0 = gex[0][0][0][0];
    break;
  default:
    std::cerr << "hoibc_ibc0::get_coeff_no_suc: mode = " << this->mode << " unknown" << std::endl;
    break;
  }
}

big_matrix<complex> hoibc::hoibc_ibc0::get_impedance(const real& k0, const std::vector<real>& f1, const std::vector<real>& f2){
  big_matrix<real> I;
  get_matrices_I(f1.size(),f2.size(),I=I);
  big_matrix<complex> impedance_ap = this->coeff.a0*I;
  return impedance_ap;
}

void hoibc::hoibc_ibc0::array_to_coeff(const std::vector<real>& x){
  this->coeff.a0 = std::complex<real>(x[0],x[1]);
}

void hoibc::hoibc_ibc0::coeff_to_array(std::vector<real>& x){
  x.clear();
  x.push_back(std::real(this->coeff.a0));
  x.push_back(std::imag(this->coeff.a0));
}

void hoibc::hoibc_ibc0::get_suc(std::vector<real>& cle, std::vector<real>& ceq, std::vector<real>& cne, std::vector<std::string>& sle, std::vector<std::string>& seq, std::vector<std::string>& sne){
  cle.clear();
  ceq.clear();
  cne.clear();
  sle.clear();
  seq.clear();
  sne.clear();
  cle.push_back(-std::real(this->coeff.a0));
  sle.push_back("-Re(a0)");
}

void hoibc::hoibc_ibc0::disp_coeff(std::ostream& out){
  out << "# Z = a0*I" << std::endl;
  print_complex(this->coeff.a0,"    a0",out);
}