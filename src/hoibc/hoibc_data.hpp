#ifndef _H_HOIBC_DATA
#define _H_HOIBC_DATA

#include "hoibc_types.hpp"

#include <iostream>
#include <vector>
#include <string>


namespace hoibc
{
  struct main_t
  {
    real frequency { 1. };          // Frequency in GHz
    std::vector<real> s1 = {0.,1.,0.};
    std::vector<real> s2 = {0.,1.,0.};
  };

  struct material_t
  {
    std::vector<complex> epsr;       // Relative permitivity
    std::vector<complex> mur;        // Relative permeability
    std::vector<real>    thickness;  // Thickness in meter
    real                 loss { 0. }; // Artificial loss to add to epsr & mur
    complex              initial_impedance[2][2] { 0. }; // The impedance to use on the first layer. If zero, it's a PEC.
  };

  struct hoibc_t
  {
    std::vector<std::string> name;
    std::vector<std::string> label;
    std::vector<char>        type;
    std::vector<bool>        suc;
    std::vector<real> inner_radius;
    std::vector<bool> normalised;
    std::vector<char> mode;
  };

  struct optim_t;

  struct data_t
  {
    main_t      main;
    material_t  material;
    hoibc_t     hoibc;
  };


  void disp_data(const data_t& data, std::ostream& out = std::cout);

  void check_data(const data_t& data);
}
#endif