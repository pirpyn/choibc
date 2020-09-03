#include "main.hpp"
#include "../hoibc/hoibc.hpp"
#include "write_impedance_errors.hpp"

// Example of a program that call the hoibc library

int main(int argc, char* argv[]) {

  assert(argc==2);
  std::cout << "# Reading data from " << argv[1] << std::endl;
  const data_out_t data_out = read_data_from_json(argv[1]);
  const hoibc::data_t data = data_out.data_t;

  // prints the parameters to stdout
  hoibc::disp_data(data);

  // A polymorphic array of all the IBC, that will store their coefficients
  // Unknown IBC are nullptr
  std::vector<hoibc::hoibc_class*> hoibc_list;

  // Compute HOIBC coefficients
  hoibc::main(data,hoibc_list);

  std::cout << std::endl;
  std::cout << "=================================================" << std::endl;
  std::cout << "HOIBC: All IBC coefficients computed with success" << std::endl;
  std::cout << "=================================================" << std::endl;
  std::cout << std::endl;

  // Write results (impedance, coeff, ...) to screen and csv files
  write_impedance_errors(data_out, hoibc_list);

  hoibc::free_hoibc_list(hoibc_list);
  return 0;
}
