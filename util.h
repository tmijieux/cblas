#ifndef TDP_UTIL_H
#define TDP_UTIL_H

#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#ifndef min
#define min(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifdef DEQUAL
#undef DEQUAL
#endif
#define DEQUAL(X_, Y_, H_) (fabs((X_) - (Y_)) < (H_))


double *tdp_matrix_new(int m/*rows*/, int n/*columns*/);
void tdp_matrix_init(int m/*rows*/, int n/*columns*/, double *mat);
void tdp_matrix_identity(int m/*row*/, int n/*column*/,
                         double *mat, int lda/*leading dimensions*/);
void tdp_matrix_print(int m/*row*/, int n/*column*/,
                      double *mat, int lda/*leading dimensions*/,
                      FILE *outstream);

double *tdp_vector_new(int m);
void tdp_vector_rand(int m, double max, double *v);
void tdp_vector_one(int m, double value, double *v);

#endif // TDP_UTIL_H
