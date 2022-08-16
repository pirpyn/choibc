#ifndef _H_HOIBC_IO
#define _H_HOIBC_IO

#include "hoibc_types.hpp"
#include <string>
#include <ostream>

namespace hoibc
{
    void set_backend(const std::string &backend);
    void disp_cmplx(std::ostream &out, const hoibc::complex &x, const std::string &s);
    const std::string get_cmt(void);
}

#ifdef _HOIBC_HAS_CPP20
#include <format>
using string_format=std::format
#else
#warning C++20 is needed to use std::format. The custom function string_format will emulate it.

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

#endif