#ifndef TDP_DDOT_H
#define TDP_DDOT_H

#include "cblas/cblas.h"

#define DEFINE_DDOT(SYMBOL)                             \
    double SYMBOL(const int N,                          \
                  const double *X, const int incX,      \
                  const double *Y, const int incY)

#define DECLARE_DDOT(SYMBOL) DEFINE_DDOT(SYMBOL);

typedef DEFINE_DDOT((*cblas_ddot_t));

#define DDOT_UNUSED_PARAMS                      \
    do {                                        \
        (void) N;                               \
        (void) X; (void) incX;                  \
        (void) Y; (void) incY;                  \
    }while(0)

#define DDOT_PARAMS  N, X, incX, Y, incY


DECLARE_DDOT(ddot_basic_Thomas);
DECLARE_DDOT(ddot_basic_Fatima_Zahra);
DECLARE_DDOT(ddot_avx_256);
DECLARE_DDOT(ddot_avxU_256); // unaligned load

DECLARE_DDOT(ddot_avx_256_fma);

#endif // TDP_DDOT_H
