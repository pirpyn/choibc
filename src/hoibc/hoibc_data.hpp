#ifndef _H_HOIBC_DATA
#define _H_HOIBC_DATA

#include "hoibc_types.hpp"

#include <iostream>
#include <string>


namespace hoibc{
  struct main_t{
    real frequency { 1. };          // Frequency in GHz
    std::vector<real> s1 = {0.,1.,0.};
    std::vector<real> s2 = {0.,1.,0.};
  };

  struct material_t{
    std::vector<complex> epsr;        // Relative permitivity
    std::vector<complex> mur;         // Relative permeability
    std::vector<real>    thickness;   // Thickness in meter
    real                 loss { 0. }; // Artificial loss to add to epsr & mur
    matrix<complex>      initial_impedance { 0. }; // The impedance to use on the first layer. If zero, it's a PEC.
  };

  struct hoibc_t{
    std::vector<std::string> name;
    std::vector<std::string> label;
    std::vector<short>       type;
    std::vector<bool>        suc;
    std::vector<real> inner_radius;
    std::vector<bool> normalised;
    std::vector<char> mode;
  };

  struct optim_t{
    integer max_iter { 100 };     // max number of iteration
    real    tol { 1.e-6 };        // tol for suc
    real    toldx  { 1e-8 };      // if `|xn+1 - xn| < toldx` then stop
    real    grad_delta { 1e-4 };  // Step to approximation gradient of constraints and cost function
    bool    no_constraints { false }; // If true, then minimsation is unconstrained
    bool    show_iter { false };
  };

  struct data_t{
    main_t      main;
    material_t  material;
    hoibc_t     hoibc;
    optim_t     optim;
  };


  void disp_data(const data_t& data, std::ostream& out = std::cout);

  void check_data(const data_t& data);
}
#endif