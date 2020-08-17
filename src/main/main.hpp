#ifndef _H_MAIN
#define _H_MAIN

#include <iostream>
#include "../hoibc/hoibc.hpp"
#include "dump_csv.hpp"
#include <string>

struct data_out_t {
    std::string basename = "hoibc";
    bool        impedance_ibc   = false;
    bool        impedance_ex    = false;
};

void write_impedance_errors(const hoibc::data_t& data, const data_out_t& data_out, std::vector<hoibc::hoibc_class*>& hoibc_list);

#endif