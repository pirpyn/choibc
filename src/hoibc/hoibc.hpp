#ifndef _H_HOIBC
#define _H_HOIBC

#include <string>

#include "hoibc_data.hpp"
#include "hoibc_class.hpp"

namespace hoibc
{

  void main(const data_t& data, std::vector<hoibc_class*>& hoibc_list);

  void setup(const data_t& data, std::vector<hoibc_class*>& hoibc_list);
  
  enum class hoibc_names;

  hoibc_names resolve_names(const std::string name);
}
#endif