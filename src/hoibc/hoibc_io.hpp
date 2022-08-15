#ifndef _H_HOIBC_IO
#define _H_HOIBC_IO

#include "hoibc_types.hpp"
#include <string>
#include <ostream>

namespace hoibc
{
    void disp_cmplx(std::ostream &out, const hoibc::complex &x, const std::string &s);
    void set_cmt(const std::string &backend);
    const std::string get_cmt(void);
    void set_fmt_cmplx(const std::string &backend);
}

#if (__cplusplus < 202002L )
#define _HOIBC_IO_FORMAT

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