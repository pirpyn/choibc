#ifndef _H_HOIBC_CLASS
#define _H_HOIBC_CLASS

#include <ostream>
#include <string>
#include <vector>
#include "hoibc_types.hpp"
#include "hoibc_data.hpp"

namespace hoibc
{
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

        virtual void get_impedance() = 0;

        virtual void array_to_coeff() = 0;

        virtual void coeff_to_array() = 0;

        virtual void get_suc() = 0;

        virtual void disp_coeff(std::ostream& out) = 0;

        void get_reflexion();

        void get_coeff(const data_t& data, const std::vector<real>& f1, const std::vector<real>& f2);

        void print_coeff(std::ostream& out = std::cout);

        void print_suc();

        void set_fourier_variables(const data_t& data, std::vector<real>& f1, std::vector<real>& f2, std::vector<real>& s1, std::vector<real>& s2);
        void set_fourier_variables(const data_t& data, std::vector<real>& f1, std::vector<real>& f2);
  };

  void print_complex(std::ostream& out, const complex& z, const std::string& name);
}
#endif