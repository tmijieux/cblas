#include <immintrin.h> //AVX
#include <assert.h>
#include "cblas.h"
#include "dgemv.h"

static
DEFINE_DGEMV(trans_mv)
{
    (void ) incX;
    (void ) incY;
    (void ) Order;
    (void ) TransA;

    int n = N/8, nr = N%8;
    double BETA[4] ALIGNED(32) = {beta, beta, beta, beta};
    __m256d b = _mm256_load_pd(BETA);

    for (int i = 0; i < n; ++i) {
        __m256d y1, y2;
        y1 = _mm256_load_pd(Y+i*8);
        y2 = _mm256_load_pd(Y+i*8+4);

        y1 = _mm256_mul_pd(y1, b);
        y2 = _mm256_mul_pd(y2, b);

        _mm256_store_pd(Y+i*8, y1);
        _mm256_store_pd(Y+i*8+4, y2);

        Y[i*8]   += alpha * ddot_avx_256(M, X, 1, A+lda*i*8, 1);
        Y[i*8+1] += alpha * ddot_avx_256(M, X, 1, A+lda*(i*8+1), 1);
        Y[i*8+2] += alpha * ddot_avx_256(M, X, 1, A+lda*(i*8+2), 1);
        Y[i*8+3] += alpha * ddot_avx_256(M, X, 1, A+lda*(i*8+3), 1);
        Y[i*8+4] += alpha * ddot_avx_256(M, X, 1, A+lda*(i*8+4), 1);
        Y[i*8+5] += alpha * ddot_avx_256(M, X, 1, A+lda*(i*8+5), 1);
        Y[i*8+6] += alpha * ddot_avx_256(M, X, 1, A+lda*(i*8+6), 1);
        Y[i*8+7] += alpha * ddot_avx_256(M, X, 1, A+lda*(i*8+7), 1);
    }

    for (int i = N-nr; i < N; ++i) {
        Y[i] = Y[i] * beta;
        Y[i] += alpha * ddot_avx_256(M, X, 1, A+lda*i, 1);
    }
}

static
DEFINE_DGEMV(notrans_mv)
{
    DGEMV_UNUSED_PARAMS;
}

DEFINE_DGEMV(dgemv_basic)
{
    assert( Order == CblasColMajor );
    assert( incX == 1 );
    assert( incY == 1 );

    if (TransA == CblasTrans)
        trans_mv(DGEMV_PARAMS_LIST);
    else
        notrans_mv(DGEMV_PARAMS_LIST);
}
