#ifndef TDP_DGEMV_H
#define TDP_DGEMV_H

#include "cblas/cblas.h"

#define DEFINE_DGEMV(SYMBOL)                                            \
    void SYMBOL(const enum CBLAS_ORDER Order,                           \
                const enum CBLAS_TRANSPOSE TransA,                      \
                const int M, const int N,                               \
                const double alpha, const double *A, const int lda,     \
                const double *X, const int incX, const double beta,     \
                double *Y, const int incY)                              \

#define DECLARE_DGEMV(SYMBOL) DEFINE_DGEMV(SYMBOL)

#define DGEMV_PARAMS_LIST                                       \
    Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY  \

#define DGEMV_UNUSED_PARAMS                     \
    do {                                        \
        (void) Order;                           \
        (void) TransA;                          \
        (void) M;                               \
        (void) N;                               \
        (void) alpha;                           \
        (void) A;                               \
        (void) lda;                             \
        (void) X;                               \
        (void) incX;                            \
        (void) beta;                            \
        (void) Y;                               \
        (void) incY;                            \
    }while(0)

typedef DEFINE_DGEMV((*cblas_dgemv_t));

DECLARE_DGEMV(dgemv_basic);

#endif // TDP_DGEMV_H
