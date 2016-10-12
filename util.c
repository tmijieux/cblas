#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "util.h"

/**
 * Pretty print the matrix 'mat'
 */
void tdp_matrix_print(int m/*rows*/, int n/*columns*/,
                      double *mat, int lda/*leading dimension*/,
                      FILE *out)
{

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j)
            fprintf(out, "%g ", mat[j*lda+i]);
        fprintf(out, "\n");
    }
    fprintf(out, "\n");
}

/**
 * Return a new zero'd (m x n) matrix
 */
double *tdp_matrix_new(int m/*rows*/, int n/*columns*/)
{
    double d;
    return calloc(m*n, sizeof d);
}

/**
 * Set the matrix elements to random values taken in interval [min, max]
 */
void tdp_matrix_rand(int m/*rows*/, int n/*columns*/,
                     double *mat, double min, double max)
{
    for (int i = 0 ; i < m*n; ++i)
        mat[i] = min + ((double)rand() / RAND_MAX)*(max-min);
}

/**
 * Set matrix elements to 0.0
 */
void tdp_matrix_zero(int m/*rows*/, int n/*columns*/, double *mat)
{
    memset(mat, 0, sizeof*mat*m*n);
}

/**
 * Set the matrix' main diagonal elements to 'value', and other to 0.0
 */
void tdp_matrix_one(int m/*rows*/, int n/*columns*/,
                    double value, double *mat, int lda/*leading dimension*/)
{
    tdp_matrix_zero(m, n, mat);
    int M = min(m, n);
    for (int j = 0; j < M; ++j)
        mat[j*lda+j] = value;
}

/**
 * Set the matrix' main diagonal elements to 'value', and other to 0.0
 */
void tdp_matrix_3one(int m/*rows*/, int n/*columns*/,
                     double v1, double v2,
                     double *mat, int lda/*leading dimension*/)
{
    tdp_matrix_zero(m, n, mat);
    int M = min(m, n);


    mat[0] = v1;
    mat[1] = v2;
    for (int j = 1; j < M-1; ++j) {
        mat[j*lda+j-1] = v2
        mat[j*lda+j] = v1;
        mat[j*lda+j+1] = v2;
    }
    mat[(M-1)*lda+(M-2)] = v2;
    mat[(M-1)*lda+(M-1)] = v1;
}


/**
 * Return new zero'd vector
 */
double *tdp_vector_new(int m)
{
    double d;
    return calloc(m, sizeof d);
}

/**
 * Set vector elements with random values taken in interval [min, max]
 */
void tdp_vector_rand(int m, double min, double max, double *v)
{
    for (int i = 0; i < m; ++i)
        v[i] = min + ((double)rand() / RAND_MAX) * (max-min);
}

/**
 * Set all vector elements to 'value'
 */
void tdp_vector_one(int m, double value, double *v)
{
    for (int i = 0; i < m; ++i)
        v[i] = value;
}

/**
 * Set vector elements to 0.0
 */
void tdp_vector_zero(int m, double *v)
{
    memset(v, 0, m*sizeof v[0]);
}

/**
 * Pretty print the vector
 */
void tdp_vector_print(int m, double *v, FILE *out)
{
    for (int i = 0; i < m; ++i)
        fprintf(out, "%g\n", v[i]);
}

#define CACHE_SIZE 25600000
void tdp_cache_garbage(void)
{
    uint64_t S = CACHE_SIZE*2;
    uint64_t s = S;
    double *a = malloc(S);
    while (s > 0) {
        int i = rand() % (S/sizeof *a);
        int k = rand() % (S/sizeof *a);
        a[i] = a[k];
        s -= sizeof *a;
    }
    free(a);
}
