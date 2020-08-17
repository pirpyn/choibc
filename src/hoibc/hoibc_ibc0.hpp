#ifndef _H_HOIBC_IBC0
#define _H_HOIBC_IBC0

#include "hoibc_class.hpp"
#include "hoibc_types.hpp"
#include <vector>

namespace hoibc
{
  
  class hoibc_ibc0 : public hoibc_class {
    public:

      struct coeff_t{
         complex a0 = {complex(0.,0.)};
      };

      struct coeff_t coeff;

      hoibc_ibc0(){};

      ~hoibc_ibc0(){};

      void get_coeff_no_suc(const std::vector<real>& f1, const std::vector<real>& f2, const big_matrix<complex>& gex, const real& k0);

      big_matrix<complex> get_impedance(const real& k0, const std::vector<real>& f1, const std::vector<real>& f2);

      void array_to_coeff(const alglib::real_1d_array& x);

      void coeff_to_array(alglib::real_1d_array& x);

      void get_suc(std::vector<real>& cle = empty_vector_real, std::vector<real>& ceq = empty_vector_real, std::vector<real>& cne = empty_vector_real, std::vector<std::string>& sle = empty_vector_string, std::vector<std::string>& seq = empty_vector_string, std::vector<std::string>& sne = empty_vector_string);
      
      void disp_coeff(std::ostream& out=std::cout);

  };

}
#endif
