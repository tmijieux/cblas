#include <stdlib.h>
#include <string.h>
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


/* Calculates the Submatrix */
void tdp_submatrix(double* A, int lda, int line, int column, int m, int n, double* B){
  for( int i=0; i<m; i++)
    for(int j=0; j<n; j++)
      B[i+lda*j]=A[i+line+lda*(j+column)];

}
