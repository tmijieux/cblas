#ifndef TDP_UTIL_H
#define TDP_UTIL_H

#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#ifndef min
#define min(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef max
#define max(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifdef DEQUAL
#undef DEQUAL
#endif
#define DEQUAL(X_, Y_, H_) (fabs((X_) - (Y_)) < (H_))

#ifdef SQUARE
#undef SQUARE
#endif
#define SQUARE(X_) ((X_)*(X_))

#ifdef CUBE
#undef CUBE
#endif
#define CUBE(X_) ((X_)*(X_)*(X_))


double *tdp_matrix_new(int m/*rows*/, int n/*columns*/);
void tdp_matrix_zero(int m/*rows*/, int n/*columns*/, double *mat);
void tdp_matrix_one(int m/*row*/, int n/*column*/,
                    double value, double *mat, int lda/*leading dimension*/);
void tdp_matrix_print(int m/*row*/, int n/*column*/,
                      double *mat, int lda/*leading dimension*/,
                      FILE *outstream);
void tdp_matrix_rand(int m/*rows*/, int n/*columns*/,
                     double *mat, double min, double max);

double *tdp_vector_new(int m);
void tdp_vector_rand(int m, double min, double max, double *v);
void tdp_vector_one(int m, double value, double *v);
void tdp_vector_zero(int m, double *v);
void tdp_vector_print(int m, double *v, FILE *out);

void tdp_print_cache_size(void);
uint64_t tdp_get_cache_size(int id);

void tdp_cache_garbage(void);

#endif // TDP_UTIL_H
