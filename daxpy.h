#ifndef TDP_DAXPY_H
#define TDP_DAXPY_H
#include "cblas/cblas.h"

#define DEFINE_DAXPY(SYMBOL)                                            \
    void SYMBOL(const int N, const double alpha, const double *X,       \
                const int incX, double *Y, const int incY)

#define DECLARE_DAXPY(SYMBOL) DEFINE_DAXPY(SYMBOL);

typedef DEFINE_DAXPY((*cblas_daxpy_t));

DECLARE_DAXPY(daxpy_basic)
DECLARE_DAXPY(daxpy_avx256);

#endif // TDP_DAXPY_H
