#ifndef _H_HOIBC
#define _H_HOIBC

#include <string>

#include "hoibc_data.hpp"
#include "hoibc_class.hpp"
#include "hoibc_constants.hpp"
#include "hoibc_math.hpp"
#include "hoibc_math_plane.hpp"

namespace hoibc
{

  enum class hoibc_names;

  void main(const data_t& data, std::vector<hoibc_class*>& hoibc_list);

  void setup(const data_t& data, std::vector<hoibc_class*>& hoibc_list);
  
  hoibc_names resolve_names(const std::string name);

  void free_hoibc_list(std::vector<hoibc_class*>& hoibc_list);
}
#endif