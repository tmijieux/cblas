#ifndef TDP_DGEMV_H
#define TDP_DGEMV_H

#include "cblas/cblas.h"

#define DEFINE_DGEMV(SYMBOL)                                            \
    void SYMBOL(const enum CBLAS_ORDER order,                           \
                const enum CBLAS_TRANSPOSE TransA,                      \
                const int M, const int N,                               \
                const double alpha,                                     \
                const double *A, const int lda,                         \
                const double *X, const int incX,                        \
                const double beta,                                      \
                double *Y, const int incY)

#define DECLARE_DGEMV(SYMBOL) DEFINE_DGEMV(SYMBOL);

typedef DEFINE_DGEMV((*cblas_dgemv_t));

DECLARE_DGEMV(dgemv_avx256);

#endif // TDP_DGEMV_H
