#include "cblas/cblas.h"

double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY)
{
    double s = 0.0;
    int i, x, y;
    for (i = x = y = 0; i != N; ++i, x += incX, y += incY)
        s += X[x] * Y[y];
    return s;
}
