#ifndef _H_HOIBC_IBC3
#define _H_HOIBC_IBC3

#include "hoibc_class.hpp"
#include "hoibc_types.hpp"
#include <vector>

namespace hoibc
{
  
  class hoibc_ibc3 : public hoibc_class {
    public:

      struct coeff_t{
         complex a0 = {complex(0.,0.)};
         complex a1 = {complex(0.,0.)};
         complex a2 = {complex(0.,0.)};
         complex b1 = {complex(0.,0.)};
         complex b2 = {complex(0.,0.)};
      };

      struct coeff_t coeff;

      hoibc_ibc3(){};

      ~hoibc_ibc3(){};

      void get_coeff_no_suc(const std::vector<real>& f1, const std::vector<real>& f2, const big_matrix<complex>& gex, const real& k0);

      big_matrix<complex> get_impedance(const real& k0, const std::vector<real>& f1, const std::vector<real>& f2);

      void array_to_coeff(const std::vector<real>& x);

      void coeff_to_array(std::vector<real>& x);

      void get_suc(std::vector<real>& cle, std::vector<real>& ceq, std::vector<real>& cne, std::vector<std::string>& sle, std::vector<std::string>& seq, std::vector<std::string>& sne);
      
      void disp_coeff(std::ostream& out=std::cout);

  };

}
#endif
