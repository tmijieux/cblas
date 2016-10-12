#include "tdp-cblas.h"
#include "tdp-ddot.h"

DEFINE_DDOT(ddot_basic_Thomas)
{
    double s = 0.0;
    int i, x, y;
    for (i = x = y = 0; i != N; ++i, x += incX, y += incY)
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
