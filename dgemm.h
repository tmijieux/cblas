#ifndef TDP_GEMM_H
#define TDP_GEMM_H

#include "cblas/cblas.h"

void cblas_dgemm_scalar(const enum CBLAS_ORDER Order,
                        const enum CBLAS_TRANSPOSE TransA,
                        const enum CBLAS_TRANSPOSE TransB,
                        const int M, const int N, const int K,
                        const double alpha, const double *A, const int lda,
                        const double *B, const int ldb, const double beta,
                        double *C, const int ldc);

#endif // TDP_GEMM_H
