#include "hoibc_ibc0.hpp"
#include "hoibc_math.hpp"
#include "hoibc_math_plane.hpp"

using namespace hoibc;

void hoibc::hoibc_ibc0::get_coeff_no_suc(const std::vector<real>& f1, const std::vector<real>& f2, const big_matrix<complex>& gex, const real& k0){
  // The IBC0 is such that a0 = Z11(0,0) = Z22(0,0)

  // We're in the vaccum
  const complex k = k0;
  const complex etar = 1.;
  const real z = 0.;

  big_matrix<complex> mAE,mBE,mAH,mBH;

  switch (this->mode){
  case 1: // if reflexion
      switch (this->type){
        case ('P'):    // For the plane we fit with respect to the reflexion coefficent
          // we only need the first element of the big matrices, so pass only first element of f1 & f2
          plane::get_matrices_AB_EH(std::vector<real>({f1[0]}),std::vector<real>({f2[0]}),k,etar,z,mAE,mBE,mAH,mBH);
          break;
        default:
          std::cerr << "hoibc_ibc0::get_coeff_no_suc: mode = " << this->mode << " type = " << this->type << "unknown" << std::endl;
          break;
      }
      this->coeff.a0 = ((mAH[0][0] + mBH[0][0]*gex[0][0]) % (mAE[0][0] + mBE[0][0]*gex[0][0]))[0][0];
        break;
  case 2: // if impedance
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

void hoibc::hoibc_ibc0::array_to_coeff(const alglib::real_1d_array& x){
  this->coeff.a0 = complex(x[0],x[1]);
}

void hoibc::hoibc_ibc0::coeff_to_array(alglib::real_1d_array& x){
  x.setlength(2);
  x[0] = std::real(this->coeff.a0);
  x[1] = std::imag(this->coeff.a0);
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
  print_complex(this->coeff.a0,"  a0",out);
}