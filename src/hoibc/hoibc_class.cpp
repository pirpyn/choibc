#include "hoibc_class.hpp"
#include "hoibc_constants.hpp"
#include "hoibc_math.hpp"
#include "hoibc_math_plane.hpp"
#include "hoibc_math_cylinder.hpp"
#include "hoibc_math_sphere.hpp"
#include <iostream>
#include <algorithm>
#include "../alglib/optimization.h"

using namespace hoibc;

big_matrix<complex> hoibc_class::get_reflexion(const real& k0, const array<real>& f1, const array<real>& f2){
  big_matrix<complex> impedance = this->get_impedance(k0,f1,f2);

  // We're in the vaccum
  // const complex k = complex(k0,0.);
  // const complex etar = complex(1.,0.);

  big_matrix<complex> reflexion;
  switch (this->type){
  case type_t::P:
    // f1 = kx, f2 = ky
    reflexion = plane::reflexion_from_impedance(f1,f2,k0,impedance);
    break;
  case type_t::C:
    // f1 = n, f2 = kz
    reflexion = cylinder::reflexion_from_impedance(f1,f2,k0,impedance,this->outer_radius);
    break;
  case type_t::S:
    // f1 = m, f2 = n
    reflexion = sphere::reflexion_from_impedance(f2,k0,impedance,this->outer_radius);
    break;
  }
  return reflexion;
}

struct costf_data_t {
  hoibc_class* ibc = NULL;
  real k0 = 0;
  std::size_t neq = 0;
  std::size_t nle = 0;
  const array<real>* f1 = NULL;
  const array<real>* f2 = NULL;
  const big_matrix<complex>* gex = NULL;
  bool no_constraints = false;
  bool show_iter = false;
};

void costf(const alglib::real_1d_array &x, alglib::real_1d_array &fi, void *ptr);
void report_iteration(const alglib::real_1d_array &x, double func, void *ptr = NULL, bool reset = false);

void hoibc_class::get_coeff(const data_t& data, const array<real>& f1, const array<real>& f2){

  // gex contains all the exact impedance or reflexion matrices for every (f1,f2) couple
  // i.e. for a plane type and mode 2 so gex[i][j][k][l] represents Z(kx[i],ky[j])_kl
  big_matrix<complex> gex;

  const real k0 = free_space_wavenumber(data.main.frequency);

  switch (this->mode){
    case mode_t::R :
      switch (this->type) {
        case type_t::P:
          gex = plane::reflexion_infinite(f1,f2,k0,data.material);
          break;
        case type_t::C :
          gex = cylinder::reflexion_infinite(f1,f2,k0,data.material,this->inner_radius);
          break;
        case type_t::S :
          gex = sphere::reflexion_infinite(f1,k0,data.material,this->inner_radius);
          break;
      }
      break;
    case mode_t::Z :
      switch (this->type) {
        case type_t::P :
          gex = plane::impedance_infinite(f1,f2,k0,data.material);
          break;
        case type_t::C :
          gex = cylinder::impedance_infinite(f1,f2,k0,data.material,this->inner_radius);
          break;
        case type_t::S :
          gex = sphere::impedance_infinite(f2,k0,data.material,this->inner_radius);
          break;
      }
      break;
  }
  this->get_coeff_no_suc(f1,f2,gex,k0);

  if (this->suc){
    std::cout << "hoibc: IBC " << this->label << ": executing the constrained optimisation algorithm ... " << std::endl;
    alglib::real_1d_array x;
    alglib::real_1d_array scale;

    // Get number of unknows from IBC, and initialise ALGLIB quantities.
    this->coeff_to_array(x);
    this->coeff_to_array(scale);

    for (alglib::ae_int_t i = 0; i < x.length(); i++){
      switch (data.optim.starting_point){
        case start_pt_t::feasible:
          // An always feasible starting point
          x[i] = 0.;
          break;
        case start_pt_t::best:
          // Already done by getting the coefficient without SUC
          break;
      }

      // We set the scale to roughly twice the no constraints solution.
      // The scale can't have zero values.
      scale[i] = abs(scale[i]) > 0. ? abs(scale[i])*2.0 : 1.0 ;
    }
    alglib::minnlcstate state;
    alglib::minnlccreatef(x, data.optim.grad_delta, state);

    // Set solver stopping criterion (epsx, maxit)
    alglib::minnlcsetcond(state, data.optim.toldx, data.optim.max_iter);

    // Set scale of unknowns
    alglib::minnlcsetscale(state, scale);

    // Use SQP solver.
    alglib::minnlcsetalgosqp(state);

    // Get number of equality and inequality constraints
    array<real> ceq, cle;
    this->get_suc(cle = cle, ceq = ceq);
    alglib::minnlcsetnlc(state, ceq.size(), cle.size());

    // Set up additionnal parameter for cost function
    costf_data_t costf_data;
    costf_data.ibc = this;
    costf_data.k0 = k0;
    costf_data.f1 = &f1;
    costf_data.f2 = &f2;
    costf_data.gex = &gex;
    costf_data.neq = ceq.size();
    costf_data.nle = cle.size();
    costf_data.no_constraints = data.optim.no_constraints;
    costf_data.show_iter = data.optim.show_iter;

    if (costf_data.no_constraints){
      std::cout << "hoibc: IBC " << this->label << " has no constraints." << std::endl;
    }
  
    // Minimize the cost function
    if (data.optim.show_iter){
      report_iteration(x,0.,NULL,true);
    }
    alglib::minnlcoptimize(state, costf, NULL, &(costf_data));

    // Get solution and termination status
    alglib::minnlcreport rep;
    alglib::minnlcresults(state, x, rep);

    this->array_to_coeff(x);
  }
}

void report_iteration(const alglib::real_1d_array &x, double func, void *ptr, bool reset){
  static int iter = 0;
  static double previous_func = func;
  static alglib::real_1d_array previous_x = x;

  if (reset){
    iter = 0;
  }

  if (iter == 0){
    // Header of the optimisation statistics
    #ifdef _HOIBC_IO_FORMAT
    std::cout << get_cmt() << " ";
    const std::string fmt_head = "%9s, %10s, %10s, ";
    std::cout << string_format(fmt_head,"iteration","f","df");
    #else
    const std::string = "{9s}, {10s}, {10s}";
    std::cout << std::format(fmt_head,"iter","f","df");
    #endif
    for (alglib::ae_int_t i = 0; i < x.length(); i++){
      #ifdef _HOIBC_IO_FORMAT
      std::cout << string_format("%10s, ",string_format("%s[%d]","x",i));
      #else
      std::cout << std::format("{10s}, ",std::format("{s}[{d}]","x",i));
      #endif
    }
    #ifdef _HOIBC_IO_FORMAT
    std::cout << string_format("%10s, ", "|dx|");
    #else
    std::cout << std::format("{10s}, ", "|dx|");
    #endif
    std::cout << std::endl;
  }
  #ifdef _HOIBC_IO_FORMAT
  const std::string fmt_head = "%9d, %10.3E, %10.3E, ";
  std::cout << string_format(fmt_head,iter,func,func-previous_func);
  #else
  const std::string = "{9d}, {10.3E}, {10.3E}";
  std::cout << std::format(fmt_head,iter,func,func-previous_func);
  #endif
  for (alglib::ae_int_t i = 0; i < x.length(); i++){
    #ifdef _HOIBC_IO_FORMAT
    std::cout << string_format("%10E, ",x[i]);
    #else
    std::cout << std::format("{10.3E}, ",x[i]);
    #endif
  }

  alglib::real_1d_array dx = x;
  alglib::vsub(dx.getcontent(),previous_x.getcontent(),dx.length());
  real norm = std::sqrt(alglib::vdotproduct(dx.getcontent(),dx.getcontent(),dx.length()));

  #ifdef _HOIBC_IO_FORMAT
  std::cout << string_format("%10.3E, ", norm );
  #else
  std::cout << std::format("{10.3E}, ", norm);
  #endif

  std::cout << std::endl;

  iter++;
  previous_func = func;
  previous_x = x;
}

void costf(const alglib::real_1d_array &x, alglib::real_1d_array &fi, void *ptr){
  costf_data_t data = *((costf_data_t*) ptr);
  data.ibc->array_to_coeff(x);

  big_matrix<complex> ap;

  switch (data.ibc->mode){
  case hoibc::mode_t::R :
    ap = data.ibc->get_reflexion(data.k0, *(data.f1), *(data.f2));
    break;
  case hoibc::mode_t::Z :
    ap = data.ibc->get_impedance(data.k0, *(data.f1), *(data.f2));
    break;
  }

  // Cost function value at current point
  fi[0] = std::pow(norm(ap - *(data.gex)),2) / std::pow(norm(*(data.gex)),2);


  if (data.show_iter){
    report_iteration(x,fi[0],NULL,false);
  }

  // Get constraints
  if (data.no_constraints) {
    for (std::size_t i = 0; i < data.neq; i++){
      fi[i+1] = 0.;
    }
    for (std::size_t i = 0; i < data.nle; i++){
      fi[i+1+data.neq] = 0.;
    }
  } else {
    array<real> cle, ceq;
    data.ibc->get_suc(cle = cle,ceq = ceq);
    for (std::size_t i = 0; i < data.neq; i++){
      fi[i+1] = ceq[i];
    }
    for (std::size_t i = 0; i < data.nle; i++){
      fi[i+1+data.neq] = cle[i];
    }
  }
}

void hoibc_class::print_coeff(std::ostream& out){
  std::noshowpos(out);
  out << "# IBC " << this->name << " type " << type_to_char(this->type) << " suc " << (this->suc ? "T": "F") << " mode " << mode_to_int(this->mode) << " (" << this->label << ")" << std::endl;
  switch (this->type) {
    case type_t::C:
    case type_t::S:
      out << get_cmt() << " inner radius " << this->inner_radius << "m, outer_radius " << this->outer_radius << "m" << std::endl;
      break;
    case type_t::P:
      break;
  }

  std::showpos(out);
  out << get_cmt();
  this->disp_coeff(out);
}

#define any(vector,logical) \
std::any_of(begin(vector), end(vector), [&tol](real x){return logical;})

#include <iomanip>      // std::setw

void hoibc_class::print_suc(const real& tol, std::ostream& out){

  array<real> cle;
  array<real> ceq;
  array<real> cne;
  array<std::string> sle;
  array<std::string> seq;
  array<std::string> sne;

  const char width = 13;

  this->get_suc(cle,ceq,cne,sle,seq,sne);

  out << "# Verifying the Sufficient Uniqueness Conditions (SUC)" << std::endl;
  out.precision(3);
  out.flags(std::ios::right | std::ios::scientific | std::ios::uppercase);
  out << std::noshowpos;
  if (cle.size()){
    if (any(cle, x <= tol)){
      out << get_cmt() << "    [OK] SUC, Negative inequality constraints, IN <= " << std::setw(10) << tol << std::endl;
      out.precision(6);
      for (std::size_t i = 0; i < cle.size(); i++){
        if (cle[i] <= tol){
          out << get_cmt() << "    IN(" << i+1 << ") = " << std::setw(width) << cle[i] << " ! " << sle[i] << std::endl;
        }
      }
    }
    if (any(cle, x > tol)){
      out << get_cmt() << "  [FAIL] SUC, Positive inequality constraints, IN > " << std::setw(10) << tol << std::endl;
      out.precision(6);
      for (std::size_t i = 0; i < cle.size(); i++){
        if (cle[i] > tol){
          out << get_cmt() << "    IN(" << i+1 << ") = " << std::setw(width) << cle[i] << " ! " << sle[i] << std::endl;
        }
      }
    }
  }

  if (ceq.size()){
    if (any(ceq, std::abs(x) <= tol)){
      out << get_cmt() << "    [OK] SUC, Zero equality constraints, |EQ| <= " << std::setw(10) << tol << std::endl;
      out.precision(6);
      for (std::size_t i = 0; i < ceq.size(); i++){
        if (ceq[i] <= tol){
          out << get_cmt() << "    EQ(" << i+1 << ") = " << std::setw(width) << ceq[i] << " ! " << seq[i] << std::endl;
        }
      }
    }
    if (any(ceq, std::abs(x) > tol)){
      out << get_cmt() << "  [FAIL] SUC, Non-zero equality constraints, |EQ| > " << std::setw(10) << tol << std::endl;
      out.precision(6);
      for (std::size_t i = 0; i+1 < ceq.size(); i++){
        if (cle[i] > tol){
          out << get_cmt() << "    EQ(" << i+1 << ") = " << std::setw(width) << ceq[i] << " ! " << seq[i] << std::endl;
        }
      }
    }
  }

  if (cne.size()){
    if (any(cne, std::abs(x) <= tol)){
      out << get_cmt() << "    [OK] SUC, Non-zero equality constraints, |NE| => " << std::setw(10) << tol << std::endl;
      out.precision(6);
      for (std::size_t i = 0; i < cne.size(); i++){
        if (cne[i] <= tol){
          out << get_cmt() << "    NE(" << i+1 << ") = " << std::setw(width) << cne[i] << " ! " << sne[i] << std::endl;
        }
      }
    }
    if (any(cne, std::abs(x) > tol)){
      out << get_cmt() << "  [FAIL] [FAIL] SUC, Zero equality constraints, |NE| < " << std::setw(10) << tol << std::endl;
      out.precision(6);
      for (std::size_t i = 0; i < cne.size(); i++){
        if (cne[i] > tol){
          out << get_cmt() << "    NE(" << i+1 << ") = " << std::setw(width) << cne[i] << " ! " << sne[i] << std::endl;
        }
      }
    }
  }
}

// Reads a data struct and creates two array with the values of the Fourier variables
void hoibc_class::set_fourier_variables(const data_t& data, array<real>& f1, array<real>& f2, array<real>& s1, array<real>& s2){
  // data: The struct to contains physics parameters,
  // f1: The 1st Fourier variable: k_x, n,   m if plane, cylinder, sphere
  // f2: The 2nd Fourier variable: k_y, k_z, n if plane, cylinder, sphere
  // s1, = f1 / free_space_wavenumber
  // s2  = f2 / free_space_wavenumber

  const real k0 = free_space_wavenumber(data.main.frequency);
  integer n1, n2;

  // This function will be called several times, so we compute at compile time the square root of two
  // for the spheric 'S' case

  switch (this->type) {
    case type_t::P :
      // f1 = kx, f2 = ky
      if (data.main.s1[2] != 0){
        n1 = static_cast<integer>((data.main.s1[1]-data.main.s1[0])/data.main.s1[2]) + 1;
      } else {
        n1 = 1;
      }
      s1 = linspace(data.main.s1[0],data.main.s1[1],n1);
      f1 = s1 * k0;

      if (data.main.s2[2] != 0){
        n2 = static_cast<integer>((data.main.s2[1]-data.main.s2[0])/data.main.s2[2]) + 1;
      } else {
        n2 = 1;
      }
      s2 = linspace(data.main.s2[0],data.main.s2[1],n2);
      f2 = s2 * k0;
      break;
    case type_t::C :
      // f1 = n, f2 = kz

      // We truncate the number of Fourier coefficient to s2*k0*outer_radius + 1
      // It should be noted that max(s2) should be at least superior or equal to 1 because Fourier coefficient in
      // front of the bessel functions Jn(k0*outer_radius) decrease exponentially as n/(k0*outer_radius) increase (same Hn)
      if (data.main.s1[1]<1){
        std::cerr << "warning: hoibc_class::set_fourier_variables: IBC " << this->label << ": enforcing max(s1) to be at least 1 to have enough Fourier coefficients ( s1 was " << data.main.s1[1] << " )" << std::endl;
      }

      n1 = static_cast<integer>(std::max(1.,data.main.s1[1])*k0*this->outer_radius) + 1;
      f1 = linspace(0,n1,n1+1);
      s1 = f1 / (k0*this->outer_radius);

      if (data.main.s2[2] != 0){
        n2 = static_cast<integer>((data.main.s2[1]-data.main.s2[0])/data.main.s2[2]) + 1;
      } else {
        n2 = 1;
      }
      s2 = linspace(data.main.s2[0],data.main.s2[1],n2);
      f2 = s2 * k0;
      break;
    case type_t::S :
      // f1 = m , f2 = n
      // The reflection doesn't depend on m, but we must provide it for the program to work
      n1 = 1;
      s1 = {0.}; // ! we set m to 0 for the correct computation of LD & LR matrices
      f1 = s1;

      // We truncate the number of Mie coefficient to s2*k0*outer_radius*sqrt(2) + 1
      // It should be noted that max(s2) should be at least superior or equal to 1 because the Mie coefficient on front
      // of the spherical bessel functions C_mn jn(k0*outer_radius) decrease exponentially as n/(k0*outer_radius) increase (same hn)
      if (data.main.s2[1] < 1.){
        std::cerr << "warning: hoibc_class::set_fourier_variables: IBC " << this->label << ": enforcing max(s2) to be at least 1 to have enough Mie coefficients ( s2 was " << data.main.s2[1] << " )" << std::endl;
      }
      n2 = static_cast<integer>(std::max(1.,data.main.s2[1])*k0*this->outer_radius*(sqrt_two)) + 1;
      f2 = linspace(0,n2,n2+1);
      s2 = f2 / (k0*this->outer_radius);
      break;
    default :
      std::cerr << "error: hoibc_class::set_fourier_variables: IBC " << this->label << " has unknown type " << type_to_char(this->type) << std::endl;
      std::exit(5);
      break;
  }
}

void hoibc_class::set_fourier_variables(const data_t& data, array<real>& f1, array<real>& f2){
  array<real> s1;
  array<real> s2;
  hoibc_class::set_fourier_variables(data,f1,f2,s1,s2);
}

void hoibc::get_matrices_I(const std::size_t& n1, const std::size_t& n2, big_matrix<real>& I, big_matrix<real>& I1, big_matrix<real>& I2){
  I.resize(0);
  I1.resize(0);
  I2.resize(0);
  matrix<real> sI;

  sI[0][0] = 1.;
  sI[1][0] = 0.;
  sI[0][1] = 0.;
  sI[1][1] = 1.;
  I = big_init<real>(n1,n2,sI);

  matrix<real> sI1;
  sI[0][0] = 1.;
  sI[1][0] = 0.;
  sI[0][1] = 0.;
  sI[1][1] = 0.;
  I1 = big_init<real>(n1,n2,sI1);

  matrix<real> sI2;
  sI[0][0] = 0.;
  sI[1][0] = 0.;
  sI[0][1] = 0.;
  sI[1][1] = 1.;
  I2 = big_init<real>(n1,n2,sI2);
}

