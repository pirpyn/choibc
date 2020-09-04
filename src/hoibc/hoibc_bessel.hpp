#ifndef _H_HOIBC_BESSEL
#define _H_HOIBC_BESSEL

#include "hoibc_types.hpp"

namespace hoibc {
    complex bessel1(const real& nu, const complex& z);
    complex bessel2(const real& nu, const complex& z);
    complex bessel1p(const real& nu, const complex& z);
    complex bessel2p(const real& nu, const complex& z);
    complex sbessel1(const real& nu, const complex& z);
    complex sbessel2(const real& nu, const complex& z);
    complex sbessel1p(const real& nu, const complex& z);
    complex sbessel2p(const real& nu, const complex& z);
}

#endif