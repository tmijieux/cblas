#include <stdlib.h>
#include <string.h>
#include "util.h"

void tdp_matrix_print(int m/*rows*/, int n/*columns*/,
                      double *mat, int lda/*leading dimensions*/,
                      FILE *outstream)
{
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j)
            fprintf(outstream, "%g ", mat[j*lda+i]);
        fprintf(outstream, "\n");
    }
    fprintf(outstream, "\n");
}

double *tdp_matrix_new(int m/*rows*/, int n/*columns*/)
{
    double d;
    return calloc(m*n, sizeof d);
}

void tdp_matrix_init(int m/*rows*/, int n/*columns*/, double *mat)
{
    memset(mat, 0, sizeof(*mat*m*n));
}

void tdp_matrix_identity(int m/*rows*/, int n/*columns*/,
                         double *mat, int lda/*leading dimension*/)
{
    tdp_matrix_init(m, n, mat);
    int M = min(m, n);
    for (int j = 0; j < M; ++j)
        mat[j*lda+j] = 1.;
}


double *tdp_vector_new(int m)
{
    double d;
    return calloc(m, sizeof d);
}

void tdp_vector_rand(int m, double max, double *v)
{
    for (int i = 0; i < m; ++i)
        v[i] = ((double)rand() / (double)RAND_MAX) * max;
}

void tdp_vector_one(int m, double value, double *v)
{
    for (int i = 0; i < m; ++i)
        v[i] = value;
}
