
#include "util.h"
#include "cblas.h"


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

    for (int i = 0; i < M; i++)
	for (int j = 0; j < N; j++) {
            C[i+j*ldc] = 0.0;
            for (int k = 0; k < K; k++)
                C[i+j*ldc] += A[k+i*lda] * B[k+j*ldb];
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
                (IS_ALIGNED(A+i*lda, 32) && IS_ALIGNED(B+j*ldb, 32))
                ? ddot_avx_256(K, A+i*lda, 1, B+j*ldb, 1)
                : ddot_avxU_256(K, A+i*lda, 1, B+j*ldb, 1);
}


DEFINE_DGEMM(dgemm_fast_sequential_beta)
{
    (void) Order;
    (void) TransA;
    (void) TransB;
    (void) alpha;
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < M; ++i)
            C[j*ldc+i] = beta * C[j*ldc+i]
                + ((IS_ALIGNED(A+i*lda, 32) && IS_ALIGNED(B+j*ldb, 32))
                   ? ddot_avx_256(K, A+i*lda, 1, B+j*ldb, 1)
                   : ddot_avxU_256(K, A+i*lda, 1, B+j*ldb, 1));
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
                (IS_ALIGNED(A+i*lda, 32) && IS_ALIGNED(B+j*ldb, 32))
                ? ddot_avx_256(K, A+i*lda, 1, B+j*ldb, 1)
                : ddot_avxU_256(K, A+i*lda, 1, B+j*ldb, 1);
}

DEFINE_DGEMM(dgemm_i)
{
    DGEMM_CHECK_PARAMS;

    for (int i = 0; i < M; i++)
	for (int j = 0; j < N; j++) {
            C[i+j*ldc] = 0.0;
            for (int k = 0; k < K; k++)
                C[i+j*ldc] += A[k+i*lda] * B[k+j*ldb];
	}
}

DEFINE_DGEMM(dgemm_j)
{
    DGEMM_CHECK_PARAMS;

    for (int j = 0; j < N; j++)
        for (int i = 0; i < M; i++){
            C[i+j*ldc] = 0.0;
            for (int k = 0; k < K; k++)
                C[i+j*ldc] += A[k+i*lda] * B[k+j*ldb];
        }
}

DEFINE_DGEMM(dgemm_k)
{
    DGEMM_CHECK_PARAMS;
    for (int j = 0; j < N; j++)
        for (int i = 0; i < M; i++)
            C[i+j*ldc] = 0.0;

    for (int k = 0; k < K; k++)
        for (int j = 0; j < N; j++)
            for (int i = 0; i < M; i++)
                C[i+j*ldc] += A[k+i*lda] * B[k+j*ldb];
}



DEFINE_DGEMM(dgemm_block)
{
    (void) alpha;
    int n, m, tb=50;

    #pragma omp parallel private(n, m) shared(tb)
    {
        n = N/tb;
        m = M/tb;

        #pragma omp for schedule(static) collapse(2)
        for (int j=0;j<m; j++) {
            for (int i=0; i<n; i++) {
                double local_beta = beta;
                for (int k = 0; k < m; k++) {
                    dgemm_fast_sequential_beta(Order,TransA,TransB,
                                               tb,tb,tb,
                                               1.0,&A[k*tb+i*tb*lda],lda,
                                               &B[k*tb+j*tb*ldb],ldb,
                                               local_beta,&C[i*tb+j*tb*ldc],ldc);
                    local_beta = 1.0;
                }
            }
        }
    } // parallel region end

    n = N%tb;
    m = M%tb;

    for (int j = 0;j < N; j++)
        for (int i = M-m; i < M; i++)
            C[j*ldc+i] = beta * C[j*ldc+i]
                + ((IS_ALIGNED(A+i*lda, 32) && IS_ALIGNED(B+j*ldb, 32))
                   ? ddot_avx_256(K, A+i*lda, 1, B+j*ldb, 1)
                   : ddot_avxU_256(K, A+i*lda, 1, B+j*ldb, 1));

    for (int j = M-m; j < M; j++)
        for (int i = 0; i < N; i++)
            C[j*ldc+i] = beta * C[j*ldc+i]
                + ((IS_ALIGNED(A+i*lda, 32) && IS_ALIGNED(B+j*ldb, 32))
                   ? ddot_avx_256(K, A+i*lda, 1, B+j*ldb, 1)
                   : ddot_avxU_256(K, A+i*lda, 1, B+j*ldb, 1));
}


