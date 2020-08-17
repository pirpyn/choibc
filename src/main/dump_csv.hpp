#ifndef _DUMP_CSH_HPP
#define _DUMP_CSH_HPP

#include <string>
#include <vector>
#include "../hoibc/hoibc.hpp"

void dump_to_csv(const std::string filename, const std::vector<hoibc::real>& f1, const std::vector<hoibc::real>& f2, const hoibc::big_matrix<hoibc::complex>& gex, const std::string& s1, const std::string& s2, const std::string& label);

// Until C++20 is available
// https://stackoverflow.com/a/26221725

#include <memory>
#include <stdexcept>

template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
    size_t size = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    std::unique_ptr<char[]> buf( new char[ size ] ); 
    snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}


#endif