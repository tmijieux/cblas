#include <assert.h>
#include "cblas/cblas.h"

void cblas_dgemm_scalar(
    const enum CBLAS_ORDER Order,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_TRANSPOSE TransB,
    const int M, const int N, const int K,
    const double alpha, const double *A, const int lda,
    const double *B, const int ldb,
    const double beta, double *C, const int ldc)
{
    assert( Order == CblasColMajor );
    assert ( alpha == 1.0 );
    assert ( beta == 0.0 );

    if (TransA == CblasTrans && TransB == CblasNoTrans) {
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < M; ++i) {
                C[j*ldc+i] = 0.0;
                for (int k = 0; k < K; ++k)
                    C[j*ldc+i] += A[i*lda+k] * B[j*ldb+k];
            }
        return;
    }
    assert( "Unsupported Transpose Configuration" && 1 == 0 );
}

//void cblas_dgemm_scalar() __attribute__((weak, alias("cblas_dgemm_scalar_mock")));


void cblas_dgemm_scalar2(const enum CBLAS_ORDER Order,
                        const enum CBLAS_TRANSPOSE TransA,
                        const enum CBLAS_TRANSPOSE TransB,
                        const int M, const int N, const int K,
                        const double alpha, const double *A, const int lda,
                        const double *B, const int ldb, const double beta,
                        double *C, const int ldc)
{
    assert( Order == CblasColMajor );
    assert ( alpha == 1.0 );
    assert ( beta == 0.0 );

    if (TransA == CblasTrans && TransB == CblasNoTrans) {
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < M; ++i)
                C[j*ldc+i] = cblas_ddot(K, A+i*lda, 1, B+j*ldb, 1);
        return;
    }
    assert( "Unsupported Transpose Configuration" && 1 == 0 );
}
