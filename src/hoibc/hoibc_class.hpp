#ifndef _H_HOIBC_CLASS
#define _H_HOIBC_CLASS

#include <ostream>
#include <string>
#include <vector>
#include "hoibc_types.hpp"
#include "hoibc_data.hpp"
#include "../alglib/ap.h"

namespace hoibc
{

  static array<real> empty_vector_real;
  static array<std::string> empty_vector_string;

  class hoibc_class {
      public:
        std::string name = "unknown";
        std::string label = "";
        type_t type = type_t::P;
        bool suc = true;
        real inner_radius = 0.;
        real outer_radius = 0.;
        bool normalised = true;
        mode_t mode = mode_t::Z;
        // struct coeff;

        inline hoibc_class(){}

        inline virtual ~hoibc_class(){}

        virtual void get_coeff_no_suc(const array<real>& f1, const array<real>& f2, const big_matrix<complex>& gex, const real& k0) = 0;

        virtual big_matrix<complex> get_impedance(const real& k0, const array<real>& f1, const array<real>& f2) = 0;

        virtual void array_to_coeff(const alglib::real_1d_array& x) = 0;

        virtual void coeff_to_array(alglib::real_1d_array& x) = 0;

        virtual void get_suc(array<real>& cle = empty_vector_real, array<real>& ceq = empty_vector_real, array<real>& cne = empty_vector_real, array<std::string>& sle = empty_vector_string, array<std::string>& seq = empty_vector_string, array<std::string>& sne = empty_vector_string) = 0;

        virtual void disp_coeff(std::ostream& out) = 0;

        big_matrix<complex> get_reflexion(const real& k0, const array<real>& f1, const array<real>& f2);

        void get_coeff(const data_t& data, const array<real>& f1, const array<real>& f2);

        void print_coeff(std::ostream& out = std::cout);

        void print_suc(const real& tol, std::ostream& out = std::cout);

        void set_fourier_variables(const data_t& data, array<real>& f1, array<real>& f2, array<real>& s1, array<real>& s2);
        void set_fourier_variables(const data_t& data, array<real>& f1, array<real>& f2);
  };

  void print_complex(const complex& z, const std::string& name, std::ostream& out = std::cout);

  void get_matrices_I(const std::size_t& n1, const std::size_t& n2, big_matrix<real>& I = empty_bigmatrix_real, big_matrix<real>& I1 = empty_bigmatrix_real, big_matrix<real>& I2 = empty_bigmatrix_real);

}
#endif