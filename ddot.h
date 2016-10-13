#ifndef TDP_DDOT_H
#define TDP_DDOT_H

#include "cblas/cblas.h"

#define DEFINE_DDOT(SYMBOL)                             \
    double SYMBOL(const int N,                          \
                  const double *X, const int incX,      \
                  const double *Y, const int incY)

#define DECLARE_DDOT(SYMBOL) DEFINE_DDOT(SYMBOL);

typedef DEFINE_DDOT((*cblas_ddot_t));


DECLARE_DDOT(ddot_basic_Thomas);
DECLARE_DDOT(ddot_basic_Fatima_Zahra);
DECLARE_DDOT(ddot_avx_256_Thomas);
DECLARE_DDOT(ddot_avx_256_fma_Thomas);

#endif // TDP_DDOT_H
