#ifndef _H_HOIBC_IBC0
#define _H_HOIBC_IBC0

#include "hoibc_class.hpp"
#include "hoibc_types.hpp"
#include <vector>
#include <map>

namespace hoibc
{
  
  class hoibc_ibc0 : public hoibc_class {
    public:

      std::map<std::string,complex> coeff{
         {"a0",{complex(0.,0.)}}
      };

      hoibc_ibc0(){};

      ~hoibc_ibc0(){};

      void get_coeff_no_suc(const std::vector<real>& f1, const std::vector<real>& f2, const big_matrix<complex>& gex, const real& k0);

      void get_impedance();

      void array_to_coeff();

      void coeff_to_array();

      void get_suc();
      
      void disp_coeff(std::ostream& out=std::cout);

  };

}
#endif
