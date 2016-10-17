#include "cblas.h"
#include "dgemv.h"
#include "util.h"

static
DEFINE_DGEMV(trans_mv)
{
    (void ) incX;
    (void ) incY;
    (void ) Order;
    (void ) TransA;

    int n = N/4, nr = N%4;
    __m256d a = _mm256_broadcast_sd(&alpha);
    __m256d b = _mm256_broadcast_sd(&beta);

    for (int i = 0; i < n; ++i) {
        __m256d d;
        double *D = (double*) &d;

        D[0] = ddot_avx_256(M, X, 1, A+lda*i*4, 1);
        D[1] = ddot_avx_256(M, X, 1, A+lda*(i*4+1), 1);
        D[2] = ddot_avx_256(M, X, 1, A+lda*(i*4+2), 1);
        D[3] = ddot_avx_256(M, X, 1, A+lda*(i*4+3), 1);

        __m256d y;
        y = _mm256_load_pd(Y+i*4);
        y = _mm256_mul_pd(y, b);
        y = MM256_FMADD_PD(d, a, y);

        _mm256_store_pd(Y+i*4, y);
    }

    for (int i = N-nr; i < N; ++i)
        Y[i] = Y[i] * beta + alpha * ddot_avx_256(M, X, 1, A+lda*i, 1);
}

static
DEFINE_DGEMV(notrans_mv)
{
    DGEMV_UNUSED_PARAMS;
    assert ( ((void)"not implemented", 1 == 0) );
}

DEFINE_DGEMV(dgemv_avx)
{
    assert( Order == CblasColMajor );
    assert( incX == 1 );
    assert( incY == 1 );

    if (TransA == CblasTrans)
        trans_mv(DGEMV_PARAMS_LIST);
    else
        notrans_mv(DGEMV_PARAMS_LIST);
}

DEFINE_DGEMV(dgemv_basic)
{
    assert( Order == CblasColMajor );
    assert( incX == 1 );
    assert( incY == 1 );
    assert( TransA == CblasTrans );

    for (int i = 0; i < N; ++i)
        Y[i] = Y[i] * beta + alpha * ddot_basic_Thomas(M, X, 1, A+lda*i, 1);
}


