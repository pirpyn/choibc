#ifndef _H_HOIBC_DATA
#define _H_HOIBC_DATA

#include "hoibc_types.hpp"

#include <iostream>
#include <string>


namespace hoibc{
  struct main_t{
    real frequency { 1. };              // Frequency in GHz
    std::array<real,3> s1 = {0.,1.,0.};  // Normalised fourier variable in the 1st dimension
    std::array<real,3> s2 = {0.,1.,0.};  // ...................................2nd dimension
  };

  struct material_t{
    array<complex> epsr;        // Relative permitivity
    array<complex> mur;         // Relative permeability
    array<real>    thickness;   // Thickness in meter
    real                 loss { 0. }; // Artificial loss to add to epsr & mur
    matrix<complex>      initial_impedance { 0. }; // The impedance to use on the first layer. If zero, it's a PEC.
  };

  // For compatibility with Fortran
  #define type_to_char(type) (type==hoibc::type_t::P?('P'):(type==hoibc::type_t::C?'C':'S'))

  enum class type_t { // Type of geometry
    P, // Plane
    C, // Cylinder
    S  // Sphere
  };

  // For compatibility with Fortran
  #define mode_to_int(mode) (mode==hoibc::mode_t::R?1:2)

  enum class mode_t { // Mode of computation
    Z, // Impedance first
    R  // Reflexion first
  };

  struct hoibc_t{
    array<std::string>  name;   // ID of the IBC
    array<std::string>  label;  // Label for output
    array<type_t>       type;   // Geometry type 'P': plane, 'C' Cylinder, 'S' Sphere
    array<bool>         suc;    // Coefficients computed w or w/o SUC
    array<real>         inner_radius; // In case of 'C' or 'S', the inner radius
    array<bool>         normalised;   // Divide L, LD, LR by k0^2
    array<mode_t>       mode;         // Computing with respect to the impedance or the reflexion
  };

  enum class start_pt_t{
    feasible,    // Starting point that verify the constraints, but may be far from the optimum
    best         // Starting point that is at an optimum, but may not verify the constraints
  };

  struct optim_t{
    integer     max_iter { 100 };                         // max number of iteration
    real        tol { 1e-6 };                             // tol for suc
    real        toldx  { 1e-8 };                          // if `|xn+1 - xn| < toldx` then stop
    real        grad_delta { 1e-4 };                      // Step to approximation gradient of constraints and cost function
    bool        no_constraints { false };                 // If true, then minimisation is unconstrained (SUC are forced to zero)
    bool        show_iter { false };                      // Display information at each iteration
    start_pt_t  starting_point { start_pt_t::feasible };  // Type of the starting point
  };

  struct data_t{
    main_t      main;
    material_t  material;
    hoibc_t     hoibc;
    optim_t     optim;
  };


  void disp_data(const data_t& data, std::ostream& out = std::cout);

  void check_data(const data_t& data);

  void check_and_set_material(integer layer_index, complex& epsr, complex&mur, complex& etar, complex& nur, const real& loss);

}
#endif