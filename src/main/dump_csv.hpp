#ifndef _DUMP_CSH_HPP
#define _DUMP_CSH_HPP

#include <string>
#include <vector>
#include "../hoibc/hoibc.hpp"

void dump_to_csv(const std::string filename, const hoibc::array<hoibc::real>& f1, const hoibc::array<hoibc::real>& f2, const hoibc::big_matrix<hoibc::complex>& gex, const std::string& s1, const std::string& s2, const std::string& label, const std::string& header);
void dump_to_csv(const std::string filename, const hoibc::array<hoibc::real>& f1, const hoibc::array<hoibc::real>& f2, const hoibc::big_matrix<hoibc::real>& gex, const std::string& s1, const std::string& s2, const std::string& label, const std::string& header);
void dump_to_csv(const std::string filename, const hoibc::array<hoibc::real>& x, const hoibc::big_matrix<hoibc::complex>& gex, const std::string& sx, const std::string& label, const std::string& header);
void dump_to_csv(const std::string filename, const hoibc::array<hoibc::real>& x, const hoibc::big_matrix<hoibc::real>& gex, const std::string& sx, const std::string& label, const std::string& header);

void set_backend(const std::string &backend);

#endif