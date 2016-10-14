//#include <zmmintrin.h> //AVX512
#include <stdio.h>
#include <immintrin.h> //AVX
#include <assert.h>

#include "cblas.h"
#include "ddot.h"

DEFINE_DDOT(ddot_basic_Thomas)
{
    double s = X[0] * Y[0];
    int i, x, y;
    for (i = 1, x = incX, y = incY; i < N; ++i, x += incX, y += incY)
        s += X[x] * Y[y];
    return s;
}

DEFINE_DDOT(ddot_avx_256)
{
    assert ( __builtin_cpu_supports("avx") );
    assert ( incX == 1 );
    assert ( incY == 1 );
    assert ( ((long)X & 31) == 0 );
    assert ( ((long)Y & 31) == 0 );

    int n = N / 4, m = N % 4;
    //printf("n=%d, m = %d\n", n, m);

    __m256d a, b, c, d;
    d = _mm256_setzero_pd();
    for (int i = 0; i < n; ++i) {
        a = _mm256_load_pd(X+i*4);
        b = _mm256_load_pd(Y+i*4);
        c = _mm256_mul_pd(a, b);
        d = _mm256_add_pd(c, d);
    }

    double s = 0.0;
    for (int i=N-m ; i < N; ++i)
    {
        s += X[i] * Y[i];
    }

    double *f = (double*)&d;
    return f[0]+f[1]+f[2]+f[3]+s;
}


DEFINE_DDOT(ddot_avxU_256)
{
    assert ( __builtin_cpu_supports("avx") );
    assert ( incX == 1 );
    assert ( incY == 1 );

    int n = N / 4, m = N % 4;
    //printf("n=%d, m = %d\n", n, m);

    __m256d a, b, c, d;
    d = _mm256_setzero_pd();
    for (int i = 0; i < n; ++i) {
        a = _mm256_loadu_pd(X+i*4);
        b = _mm256_loadu_pd(Y+i*4);
        c = _mm256_mul_pd(a, b);
        d = _mm256_add_pd(c, d);
    }

    double s = 0.0;
    for (int i=N-m ; i < N; ++i)
    {
        s += X[i] * Y[i];
    }

    double *f = (double*)&d;
    return f[0]+f[1]+f[2]+f[3]+s;
}

DEFINE_DDOT(ddot_avx_256_fma)
{
    #ifdef __FMA__
    assert ( ((long)X & 31) == 0 );
    assert ( ((long)Y & 31) == 0 );

    assert ( incX == 1 );
    assert ( incY == 1 );

    int n = N / 4, m = N % 4;

    __m256d a, b, c;
    c = _mm256_setzero_pd();
    for (int i = 0; i < n; ++i) {
        a = _mm256_load_pd(X+i*4);
        b = _mm256_load_pd(Y+i*4);
        c = _mm256_fmadd_pd(a, b, c);
    }

    double s = 0.0;
    for (int i=N-m ; i < N; ++i)
        s += X[i] * Y[i];

    double *f = (double*)&c;
    return f[0]+f[1]+f[2]+f[3]+s;
    #else
    DDOT_UNUSED_PARAMS;
    fprintf(stderr, "FMA not supported!!");
    return 0.0;
    #endif
}

DEFINE_DDOT(ddot_basic_Fatima_Zahra)
{
    double Z=0;
    for(int i=0;i<N;i++)
        Z=Z+X[i*incX]*Y[i*incY];
    return Z;
}
