#ifndef _H_READ_JSON
#define _H_READ_JSON

#include "../hoibc/hoibc.hpp"

struct data_out_t {
    hoibc::data_t data_t;
    std::string basename = "hoibc";
    bool        impedance_ibc   = false;
    bool        impedance_ex    = false;
    bool        impedance_err   = false;
    bool        coeff           = false;
    bool        reflexion_ibc   = false;
    bool        reflexion_ex    = false;
};

data_out_t read_data_from_json(const std::string& filename);

#endif