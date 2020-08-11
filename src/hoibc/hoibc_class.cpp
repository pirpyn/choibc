#include "hoibc_class.hpp"
#include "hoibc_constants.hpp"
#include "hoibc_math.hpp"
#include "hoibc_math_plane.hpp"
#include "hoibc_math_cylinder.hpp"
#include "hoibc_math_sphere.hpp"
#include <iostream>
#include <algorithm>
//#include <cstdlib>

using namespace hoibc;
using std::vector;

void hoibc_class::get_reflexion(){
  std::cout << "hoibc_class::get_reflexion: i do nothing" << std::endl;
}

void hoibc_class::get_coeff(const data_t& data, const vector<real>& f1, const vector<real>& f2){

  // gex contains all the exact impedance or reflexion matrices for every (f1,f2) couple
  // i.e. for a plane type and mode 2 so gex[i][j][k][l] represents Z(kx[i],ky[j])_kl
  big_matrix<complex> gex;

  const real k0 = free_space_impedance(data.main.frequency);

  switch (this->mode){
    case 1 :
      switch (this->type) {
        case 'P' :
          reflexion_infinite_plane();
          break;
        case 'C' :
          reflexion_infinite_cylinder();
          break;
        case 'S' :
          reflexion_infinite_sphere();
          break;
      }
      break;
    case 2 :
      switch (this->type) {
        case 'P' :
          gex = impedance_infinite_plane(f1,f2,k0,data.material);
          break;
        case 'C' :
          impedance_infinite_cylinder();
          break;
        case 'S' :
          impedance_infinite_sphere();
          break;
      }
      break;
  }
  this->get_coeff_no_suc(f1,f2,gex,k0);

  std::cout << "hoibc_class::get_coeff: not finished" << std::endl;
}

/*
    ! Then if SUC, use slsqp
    if (this->suc) then

      write(output_unit,'(a)',advance='no') '# Executing the constrained optimisation algorithm ... '

      ! Allocate x with the right number of element, depending on the HOIBC
      call this->coeff_to_array(x)
      n = size(x)
      ! Get the number of constraints
      call this->get_suc(cle,ceq,cne,sle,seq,sne)
      meq = size(ceq)
      m = meq + size(cle)


      allocate(xl(n))
      allocate(xu(n))

      ! Set bounds to big values ( heuristic )
      ! TODO Find an other constrained algorithm that doesn't rely on bounds
      xl(:) = -10._wp
      xu(:) = 10._wp

      ! Its mandatory to start inside the constraints for slsqp to converge
      ! zero is part of this space luckily, so let start from it
      x(:) = 0.

      ! The cost function to minimise.
      f => costf
      ! At the moment we evaluate only approximated gradient
      g => gradf

      if (.not.allocated(g_ex)) then
        write(error_unit,'(*(a))') 'error:get_coeff: forgot to set g_ex for type ',this->type,' ibc'
        error stop
      end if
      if (allocated(g_f1)) deallocate(g_f1)
      if (allocated(g_f2)) deallocate(g_f2)
      if (allocated(g_ibc)) deallocate(g_ibc)
 
      if (allocated(g_xlast)) deallocate(g_xlast)
      allocate(g_xlast(n))

      g_k0 = k0
      allocate(g_f1(size(f1)))
      g_f1 = f1
      allocate(g_f2(size(f2)))
      g_f2 = f2
      g_ibc = self
      g_gradient_delta = optim%grad_delta
      g_no_constraints = optim%no_constraints

      if (g_no_constraints) then
        write(output_unit,'(3(a))') 'IBC ',this->name,' has no constraints'
      end if

      if (optim%show_iter) then
        report => report_iteration
      else
        report => empty_report
      end if

      call solver%initialize(n,m,meq,optim%max_iter,optim%acc,f,g,xl,xu, &
        linesearch_mode=optim%linesearch_mode, &
        status_ok=status_ok,report=report,toldf=optim%toldf,toldx=optim%toldx)

      if (.not.status_ok) then
        error stop 'error calling slsqp%initialize.'
      end if

      call solver%optimize(x,istat,iterations)

      if (istat.ne.0) then
        write(*,'(a,1x,i0)') 'istat      :', istat
        error stop 'Error during slsqp'
      end if

      if (optim%show_iter) then
        allocate(final_c(m))
        call costf(solver,x,final_f,final_c)
        call report_iteration(solver,iterations,x,final_f,final_c)
        deallocate(final_c)
      end if

      deallocate(g_ex)
      deallocate(g_f1)
      deallocate(g_f2)
      deallocate(g_ibc)
      deallocate(xl)
      deallocate(xu)
      deallocate(g_xlast)
      ! Save the solution in the IBC
      call this->array_to_coeff(x)
    end if
  end subroutine
  */

void hoibc_class::print_coeff(std::ostream& out){
  std::noshowpos(out);
  out << "# IBC " << this->name << " type " << this->type << " suc " << (this->suc ? "T": "F") << " mode " << this->mode << " (" << this->label << ")" << std::endl;
  std::showpos(out);
  this->disp_coeff(out);
}

void hoibc_class::print_suc(){
  std::cout << "hoibc_class::get_print_suc: i do nothing" << std::endl;
}

// Reads a data struct and creates two array with the values of the Fourier variables 
void hoibc_class::set_fourier_variables(const data_t& data, vector<real>& f1, vector<real>& f2, vector<real>& s1, vector<real>& s2){
  // data: The struct to contains physics parameters, 
  // f1: The 1st Fourier variable: k_x, n,   m if plane, cylinder, sphere
  // f2: The 2nd Fourier variable: k_y, k_z, n if plane, cylinder, sphere
  // s1, = f1 / free_space_wavenumber
  // s2  = f2 / free_space_wavenumber */

  const real k0 = hoibc::free_space_impedance(data.main.frequency);
  integer n1, n2;

  // This function will be called several times, so we compute at compile time the square root of two
  // for the spheric 'S' case

  switch (this->type) {
    case 'P' :
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
    case 'C' :
      // f1 = n, f2 = kz

      // We truncate the number of Fourier coefficient to s2*k0*outer_radius + 1
      // It should be noted that max(s2) should be at least superior or equal to 1 because Fourier coefficient in 
      // front of the bessel functions Jn(k0*outer_radius) decrease exponentially as n/(k0*outer_radius) increase (same Hn)
      if (data.main.s1[1]<1) {
        std::cerr << "warning: hoibc_class::set_fourier_variables: enforcing max(s1) to be at least 1 to have enough Fourier coefficients ( s1 was " << data.main.s1[1] << " )" << std::endl;
      }
      
      n1 = static_cast<integer>(std::max(1.,data.main.s1[1])*k0*this->outer_radius) + 1;
      f1 = linspace(0,n1,n1+1);
      s1 = f1 / (k0*this->outer_radius);

      n2 = int((data.main.s2[1]-data.main.s2[0])/data.main.s2[2]) + 1;
      s2 = linspace(data.main.s2[0],data.main.s2[1],n2);
      f2 = s2 * k0;
      break;
    case 'S' :
      // f1 = m , f2 = n
      // The reflection doesn't depend on m, but we must provide it for the program to work
      n1 = 1;
      s1 = {0.}; // ! we set m to 0 for the correct computation of LD & LR matrices
      f1 = s1;

      // We truncate the number of Mie coefficient to s2*k0*outer_radius*sqrt(2) + 8 (same as Matlab)
      // It should be noted that max(s2) should be at least superior or equal to 1 because the Mie coefficient on front
      // of the spherical bessel functions jn(k0*outer_radius) decrease exponentially as n/(k0*outer_radius) increase (same hn)
      if (data.main.s2[1] < 1.) {
        std::cerr << "warning: hoibc_class::set_fourier_variables: enforcing max(s2) to be at least 1 to have enough Mie coefficients ( s2 was " << data.main.s2[1] << " )" << std::endl;
      }
      n2 = static_cast<integer>(std::max(1.,data.main.s2[1])*k0*this->outer_radius*(sqrt_two)) + 8;
      f2 = linspace(0,n2,n2+1);
      s2 = f2 / (k0*this->outer_radius);
      break;
    default :
      std::cerr << "error: hoibc_class::set_fourier_variables: IBC has unknown type " << this->type << std::endl;
      std::exit(5);
      break;
  }
}

void hoibc_class::set_fourier_variables(const data_t& data, vector<real>& f1, vector<real>& f2){
  vector<real> s1;
  vector<real> s2;
  hoibc_class::set_fourier_variables(data,f1,f2,s1,s2);
}

void hoibc::print_complex(std::ostream& out, const complex& z, const std::string& name){
  out << name << " = " << z << std::endl;
}
