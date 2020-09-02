#ifndef _H_MAIN
#define _H_MAIN

#include "../hoibc/hoibc.hpp"
#include "read_json.hpp"
#include <string>

void write_impedance_errors(const data_out_t& data_out, std::vector<hoibc::hoibc_class*>& hoibc_list);

#endif