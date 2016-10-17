#include "util.h"
#include "cblas.h"
#include "daxpy.h"

#define DAXPY_CHECK_INC_1                       \
    do {                                        \
        assert( incX == 1 );                    \
        assert( incY == 1 );                    \
    }while(0)

#define DAXPY_CHECK_ALIGNED_INPUT               \
    do {                                        \
        assert( IS_ALIGNED(X, 32) );            \
        assert( IS_ALIGNED(Y, 32) );            \
    }while(0)

DEFINE_DAXPY(daxpy_basic)
{
    for (int i = 0; i < N; ++i)
        Y[i*incY] += alpha * X[i*incX];
}

DEFINE_DAXPY(daxpy_basic1)
{
    DAXPY_CHECK_INC_1;

    for (int i = 0; i < N; ++i)
        Y[i] += alpha * X[i];
}

DEFINE_DAXPY(daxpy_avx256)
{
    DAXPY_CHECK_INC_1;
    DAXPY_CHECK_ALIGNED_INPUT;

    int n = N/4;
    int m = N%4;

    __m256d a = _mm256_broadcast_sd(&alpha);
    for (int i = 0; i < n; ++i) {
        __m256d x, y;
        x = _mm256_load_pd(X+i*4);
        y = _mm256_load_pd(Y+i*4);
        y = MM256_FMADD_PD(a, x, y);
        _mm256_store_pd(Y+i*4, y);
    }
    for (int i = N-m; i < N; ++i)
        Y[i] += alpha * X[i];
}

