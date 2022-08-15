#include "hoibc_io.hpp"

using namespace hoibc;

static std::string cmt = "";
static std::string fmt_cmplx = "";

const std::string hoibc::get_cmt(void)
{
    return cmt;
}

void hoibc::disp_cmplx(std::ostream &out, const hoibc::complex &x, const std::string &s)
{
    out << string_format(fmt_cmplx,s.c_str(),x.real(),x.imag()) << std::endl;
}

void hoibc::set_fmt_cmplx(const std::string &backend)
{
    if (backend == "csv" || backend == "comma" || backend == "paraview" || backend == "libreoffice" ){
        fmt_cmplx = "  %s,%15.8E,%15.8E";
    }
    else if (backend == "semicolon"){
        fmt_cmplx = "  %s;%15.8E;%15.8E";
    }
    else if (backend == "scilab"){
        fmt_cmplx = "  %s = %15.8E + %15.8E*%%i";
    }
    else if (backend == "octave" || backend == "matlab"){
        fmt_cmplx = "  %s = %15.8E + %15.8E*1i";
    }
    else{
        fmt_cmplx = "  %s = ( %15.8E, %15.8E )";
    }
}

void hoibc::set_cmt(const std::string &backend)
{
    if (backend == "matlab" || backend == "octave"){
        cmt = "%";
    }
    else if ( backend == "scilab" ){
        cmt = "//";
    }
    else{
        cmt = "#";
    }
}
