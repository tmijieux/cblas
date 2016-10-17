
#include "util.h"
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

DEFINE_DDOT(ddot_basic_Fatima_Zahra)
{
    double Z=0;
    for(int i=0;i<N;i++)
        Z=Z+X[i*incX]*Y[i*incY];
    return Z;
}

DEFINE_DDOT(ddot_avx_256)
{
    assert ( __builtin_cpu_supports("avx") );
    assert ( incX == 1 );
    assert ( incY == 1 );
    assert ( IS_ALIGNED(X, 32) );
    assert ( IS_ALIGNED(Y, 32) );

    int n = N / 4, m = N % 4;
    //printf("n=%d, m = %d\n", n, m);

    __m256d a, b, c;
    c = _mm256_setzero_pd();
    for (int i = 0; i < n; ++i) {
        a = _mm256_load_pd(X+i*4);
        b = _mm256_load_pd(Y+i*4);
        c = MM256_FMADD_PD(a, b, c);
    }

    double s = 0.0;
    for (int i=N-m ; i < N; ++i)
        s += X[i] * Y[i];

    double *f = (double*)&c;
    return f[0]+f[1]+f[2]+f[3]+s;
}

DEFINE_DDOT(ddot_avxU_256)
{
    assert ( __builtin_cpu_supports("avx")  );
    assert ( incX == 1 );
    assert ( incY == 1 );

    int n = N / 4, m = N % 4;

    __m256d a, b, c;
    c = _mm256_setzero_pd();
    for (int i = 0; i < n; ++i) {
        a = _mm256_loadu_pd(X+i*4);
        b = _mm256_loadu_pd(Y+i*4);
        c = MM256_FMADD_PD(a, b, c);
    }

    double s = 0.0;
    for (int i = N-m ; i < N; ++i)
        s += X[i] * Y[i];

    double *f = (double*)&c;
    return f[0]+f[1]+f[2]+f[3]+s;
}
