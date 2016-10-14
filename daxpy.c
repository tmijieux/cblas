#include <immintrin.h> //AVX
#include <assert.h>
#include "cblas.h"
#include "daxpy.h"

DEFINE_DAXPY(daxpy_basic)
{
    for (int i = 0; i < N; ++i)
        Y[i*incY] += alpha * X[i*incX];
}

#define DAXPY_CHECK_PARAMS                      \
    do {                                        \
        assert( incX == 1 );                    \
        assert( incY == 1 );                    \
    }while(0)


DEFINE_DAXPY(daxpy_avx256)
{
    DAXPY_CHECK_PARAMS;

    assert( IS_ALIGNED(X, 32) );
    assert( IS_ALIGNED(Y, 32) );

    int n = N/4;
    int m = N%4;

    __m256d a, x, y;
    double ALPHA[4] ALIGNED(32) = {alpha, alpha, alpha, alpha};
    a = _mm256_load_pd(ALPHA);

    x = _mm256_setzero_pd();

    for (int i = 0; i < n; ++i) {
        x = _mm256_loadu_pd(X+i*4);
        y = _mm256_loadu_pd(Y+i*4);
        #if __FMA__
        y = _mm256_fmadd_pd(a, x, y);
        #else
        y = _mm256_add_pd(y, _mm256_mul_pd(a, x));
        #endif
        _mm256_store_pd(Y+i*4, y);
    }
    for (int i = N-m-1; i < N; ++i)
        Y[i] += alpha * X[i];
}


