#ifndef TDP_GEMM_H
#define TDP_GEMM_H

#include "cblas/cblas.h"

#define DEFINE_DGEMM(SYMBOL)                            \
    void SYMBOL(const enum CBLAS_ORDER Order,           \
                const enum CBLAS_TRANSPOSE TransA,      \
                const enum CBLAS_TRANSPOSE TransB,      \
                const int M, const int N,const int K,   \
                const double alpha,                     \
                const double *A,const int lda,          \
                const double *B, const int ldb,         \
                const double beta,                      \
                double *C, const int ldc)

#define DECLARE_DGEMM(SYMBOL) DEFINE_DGEMM(SYMBOL);

typedef DEFINE_DGEMM((*cblas_dgemm_t));

DECLARE_DGEMM(dgemm_scalar_Fatima_Zahra);
DECLARE_DGEMM(dgemm_scalar_Thomas);
DECLARE_DGEMM(dgemm_scalar2_Thomas);
DECLARE_DGEMM(dgemm_i);
DECLARE_DGEMM(dgemm_j);
DECLARE_DGEMM(dgemm_k);
DECLARE_DGEMM(dgemm_fast_sequential);
DECLARE_DGEMM(dgemm_fast_OMP);

#endif // TDP_GEMM_H
