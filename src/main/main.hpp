#ifndef _H_MAIN
#define _H_MAIN

#include <iostream>
#include "../hoibc/hoibc.hpp"

void write_impedance_errors(const hoibc::data_t& data, std::vector<hoibc::hoibc_class*>& hoibc_list);

// Until C++20 is available
// https://stackoverflow.com/a/26221725

#include <memory>
#include <string>
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