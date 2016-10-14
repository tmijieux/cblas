#include <assert.h>

#include "dgemm.h"
#include "ddot.h"

#define DGEMM_CHECK_PARAMS                      \
    do {                                        \
        assert( TransA == CblasTrans );         \
        assert( TransB == CblasNoTrans );       \
        assert( alpha == 1.0 );                 \
        assert( beta == 0.0 );                  \
        assert( Order == CblasColMajor );       \
    }while(0);

DEFINE_DGEMM(dgemm_scalar_Fatima_Zahra)
{
    DGEMM_CHECK_PARAMS;

    for (int i=0; i<M; i++)
	for (int j=0; j<N; j++) {
            C[i+j*ldc] = 0.0;
            for (int k=0; k<K; k++)
                C[i+j*ldc]=C[i+j*ldc]+A[k+i*lda]*B[k+j*ldb];
	}
}

DEFINE_DGEMM(dgemm_scalar_Thomas)
{
    DGEMM_CHECK_PARAMS;

    for (int j = 0; j < N; ++j)
        for (int i = 0; i < M; ++i) {
            C[j*ldc+i] = 0.0;
            for (int k = 0; k < K; ++k)
                C[j*ldc+i] += A[i*lda+k] * B[j*ldb+k];
        }
}

DEFINE_DGEMM(dgemm_scalar2_Thomas)
{
    DGEMM_CHECK_PARAMS;

    for (int j = 0; j < N; ++j)
        for (int i = 0; i < M; ++i)
            C[j*ldc+i] = ddot_basic_Thomas(K, A+i*lda, 1, B+j*ldb, 1);
}

DEFINE_DGEMM(dgemm_fast_sequential)
{
    DGEMM_CHECK_PARAMS;

    for (int j = 0; j < N; ++j)
        for (int i = 0; i < M; ++i)
            C[j*ldc+i] =
                ((((long)(A+i*lda) & 31) == 0) && (((long)(B+j*ldb) & 31) == 0))
                ? ddot_avx_256(K, A+i*lda, 1, B+j*ldb, 1)
                : ddot_avxU_256(K, A+i*lda, 1, B+j*ldb, 1);
}

DEFINE_DGEMM(dgemm_OMP)
{
    DGEMM_CHECK_PARAMS;

    #pragma omp parallel for schedule(static) collapse(2)
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < M; ++i) {
            C[j*ldc+i] = 0.0;
            for (int k = 0; k < K; ++k)
                C[j*ldc+i] += A[i*lda+k] * B[j*ldb+k];
        }
}

DEFINE_DGEMM(dgemm_fast_OMP)
{
    DGEMM_CHECK_PARAMS;

    #pragma omp parallel for schedule(static) collapse(2)
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < M; ++i)
            C[j*ldc+i] =
                ((((long)(A+i*lda) & 31) == 0) && (((long)(B+j*ldb) & 31) == 0))
                ? ddot_avx_256(K, A+i*lda, 1, B+j*ldb, 1)
                : ddot_avxU_256(K, A+i*lda, 1, B+j*ldb, 1);
}

DEFINE_DGEMM(dgemm_i)
{
    DGEMM_CHECK_PARAMS;

    for (int i=0; i<M; i++)
	for (int j=0; j<N; j++) {
            C[i+j*ldc]=0;
            for (int k=0; k<K; k++)
                C[i+j*ldc]=C[i+j*ldc]+A[k+i*lda]*B[k+j*ldb];
	}
}

DEFINE_DGEMM(dgemm_j)
{
    DGEMM_CHECK_PARAMS;

    for (int j=0; j<N; j++)
        for (int i=0; i<M; i++){
            C[i+j*ldc]=0;
            for (int k=0; k<K; k++)
                C[i+j*ldc]=C[i+j*ldc]+A[k+i*lda]*B[k+j*ldb];
        }
}

DEFINE_DGEMM(dgemm_k)
{
    DGEMM_CHECK_PARAMS;
    for (int j=0; j<N; j++)
        for (int i=0; i<M; i++)
            C[i+j*ldc]=0.0;

    for (int j=0;j<N; j++)
        for (int i=0; i<M; i++)
            C[i+j*ldc] = 0.0;

    for (int k=0; k<K; k++)
        for (int j=0;j<N; j++)
            for (int i=0; i<M; i++)
                C[i+j*ldc]=C[i+j*ldc]+A[k+i*lda]*B[k+j*ldb];
}


