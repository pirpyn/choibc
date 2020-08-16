#ifndef _H_HOIBC_CLASS
#define _H_HOIBC_CLASS

#include <ostream>
#include <string>
#include <vector>
#include "hoibc_types.hpp"
#include "hoibc_data.hpp"
namespace hoibc
{

  using suc = std::vector<real,std::string>;

  class hoibc_class {
      public:
        std::string name = "unknown";
        std::string label = "";
        char type = 'P';
        bool suc = true;
        real inner_radius = 0.;
        real outer_radius = 0.;
        bool normalised = true;
        short mode = 1;
        // struct coeff;

        inline hoibc_class(){}

        inline virtual ~hoibc_class(){}

        virtual void get_coeff_no_suc(const std::vector<real>& f1, const std::vector<real>& f2, const big_matrix<complex>& gex, const real& k0) = 0;

        virtual big_matrix<complex> get_impedance(const real& k0, const std::vector<real>& f1, const std::vector<real>& f2) = 0;

        virtual void array_to_coeff(const std::vector<real>& x) = 0;

        virtual void coeff_to_array(std::vector<real>& x) = 0;

        virtual void get_suc(std::vector<real>& cle, std::vector<real>& ceq, std::vector<real>& cne, std::vector<std::string>& sle, std::vector<std::string>& seq, std::vector<std::string>& sne) = 0;

        virtual void disp_coeff(std::ostream& out) = 0;

        big_matrix<complex> get_reflexion(const real& k0, std::vector<real>& f1, std::vector<real>& f2);

        void get_coeff(const data_t& data, const std::vector<real>& f1, const std::vector<real>& f2);

        void print_coeff(std::ostream& out = std::cout);

        void print_suc(const real& tol, std::ostream& out = std::cout);

        void set_fourier_variables(const data_t& data, std::vector<real>& f1, std::vector<real>& f2, std::vector<real>& s1, std::vector<real>& s2);
        void set_fourier_variables(const data_t& data, std::vector<real>& f1, std::vector<real>& f2);
  };

  void print_complex(const complex& z, const std::string& name, std::ostream& out = std::cout);

  void get_matrices_I(const std::size_t& n1, const std::size_t& n2, big_matrix<real>& I = empty_bigmatrix_real, big_matrix<real>& I1 = empty_bigmatrix_real, big_matrix<real>& I2 = empty_bigmatrix_real);

}
#endif