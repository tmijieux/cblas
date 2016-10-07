#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <string.h>

#include "util.h"
#include "cblas/cblas.h"
#include "perf/perf.h"

void test_matrix_print(void)
{
    double M[4 * 6];
    M[0] = 11; M[4] = 12; M[8] = 13;  M[12] = 14; M[16] = 15; M[20] = 16;
    M[1] = 21; M[5] = 22; M[9] = 23;  M[13] = 24; M[17] = 25; M[21] = 26;
    M[2] = 31; M[6] = 32; M[10] = 33; M[14] = 34; M[18] = 35; M[22] = 36;
    M[3] = 41; M[7] = 42; M[11] = 43; M[15] = 44; M[19] = 45; M[23] = 46;

    tdp_matrix_print(4, 6, M, 4, stdout);
    tdp_matrix_print(2, 2, M+9, 4, stdout);
}

void test_matrix_allocate(void)
{
    double *m = tdp_matrix_new(10, 10);
    tdp_matrix_print(10, 10, m, 10, stdout);
    free(m);
}

void test_matrix_identity(void)
{
    double *m = tdp_matrix_new(10, 10);
    tdp_matrix_identity(10, 10, m, 10);
    tdp_matrix_print(10, 10, m, 10, stdout);
    free(m);

    m = tdp_matrix_new(2, 10);
    tdp_matrix_identity(2, 10, m, 2);
    tdp_matrix_print(2, 10, m, 2, stdout);
    free(m);
}

void test_vector_ddot(void)
{
    double X[] = { 1, 2, 3, 4, 5, 6 };
    assert( DEQUAL( cblas_ddot(6, X, 1, X, 1), 91.0, 0.01) );
}

#define NB_ITER 1000

void bench_vector_ddot_incONE(void)
{
    int m = 50;
    while ( m < 1000000 ) {
        double *v1, *v2;
        v1 = tdp_vector_new(m); v2 = tdp_vector_new(m);
        tdp_vector_rand(m, DBL_MAX, v1); tdp_vector_rand(m, 500.0, v2);

        perf_t p1, p2;
        perf(&p1);
        for (int i = 0; i < NB_ITER; ++i)
            cblas_ddot(m, v1, 1, v2, 1);
        perf(&p2);

        perf_diff(&p1, &p2);
        printf("m = %6d ", m);
        uint64_t nb_op = 2 * m * NB_ITER;
        printf("Mflops = %8g | time(µs) = ", perf_mflops(&p2, nb_op));
        perf_printmicro(&p2);

        free(v1); free(v2);
        m = (int) (1.25 * m);
    }
}


void test_dgemm_scalar(void)
{
    int m = 4, n = 5, k = 3;
    double A[m*k], B[k*n], C[m*n], D[m*n];

    A[0] = 1.; A[4] = 8.; A[8] = 9.;
    A[1] = 2.; A[5] = 7.; A[9] = 10.;
    A[2] = 3.; A[6] = 6.; A[10] = 11.;
    A[3] = 4.; A[7] = 5.; A[11] = 12.;

    B[0] = 11; B[3] = 12; B[6] = 13; B[9] = 14;  B[12] = 15;
    B[1] = 21; B[4] = 22; B[7] = 23; B[10] = 24; B[13] = 25;
    B[2] = 31; B[5] = 32; B[8] = 33; B[11] = 34; B[14] = 35;

    D[0] = 11; D[4] = 12; D[8] = 13; D[12] = 13; D[16] = 13;
    D[1] = 21; D[5] = 22; D[9] = 23; D[13] = 13; D[17] = 13;
    D[2] = 31; D[6] = 32; D[10]= 33; D[14] = 13; D[18] = 13;
    D[3] = 41; D[7] = 42; D[11]= 43; D[15] = 13; D[19] = 13;

    cblas_dgemm_scalar(CblasColMajor, CblasTrans, CblasNoTrans,
                       m, n, k, 0.0, A, m, B, k, 0.0, C, m);

    assert( memcmp(C, D, 4*5*sizeof C[0]) == 0 );
}

int main(int argc, char **argv)
{
    (void) argv;

    srand(time(NULL) + (long)&argc);
    test_matrix_print();
    test_matrix_allocate();
    test_matrix_identity();

    test_vector_ddot();

    test_dgemm_scalar();

    bench_vector_ddot_incONE();

    return EXIT_SUCCESS;
}
