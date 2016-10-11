#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <string.h>

#include "util.h"
#include "dgemm.h"
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

void test_matrix_one(void)
{
    double *m = tdp_matrix_new(10, 10);
    tdp_matrix_one(10, 10, 12.0, m, 10);
    tdp_matrix_print(10, 10, m, 10, stdout);
    free(m);

    m = tdp_matrix_new(2, 10);
    tdp_matrix_one(2, 10, -42.0, m, 2);
    tdp_matrix_print(2, 10, m, 2, stdout);
    free(m);

    m = tdp_matrix_new(7, 4);
    tdp_matrix_one(7, 4, -13.0, m, 7);
    tdp_matrix_print(7, 4, m, 7, stdout);
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
    printf("\n");
    int m = 50;
    while ( m < 1000000 ) {
        double *v1, *v2;
        v1 = tdp_vector_new(m); v2 = tdp_vector_new(m);
        tdp_vector_rand(m, 42., DBL_MAX, v1); tdp_vector_rand(m, -37.0, 500.0, v2);

        tdp_cache_garbage();
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
        m *= 1.25;
    }
}

void bench_dgemm_scalar(void)
{
    printf("\n");
    for (int m = 100; m <= 1000; m += 50) {
        double *M1, *M2, *M3;
        M1 = tdp_matrix_new(m, m);
        M2 = tdp_matrix_new(m, m);
        M3 = tdp_matrix_new(m, m);

        tdp_matrix_rand(m, m, M1, 42., DBL_MAX);
        tdp_matrix_rand(m, m, M2, -37.0, 500.0);

        tdp_cache_garbage();
        perf_t p1, p2;
        perf(&p1);
        cblas_dgemm_scalar(CblasColMajor, CblasTrans, CblasNoTrans,
                           m, m, m, 1.0, M1, m, M2, m, 0.0, M3, m);
        perf(&p2);

        perf_diff(&p1, &p2);
        printf("m = %6d ", m);
        uint64_t nb_op = 2 * m * m * m;
        printf("Mflops = %8g | time(µs) = ", perf_mflops(&p2, nb_op));
        perf_printmicro(&p2);

        free(M1); free(M2); free(M3);
    }
}

void test_dgemm_scalar(void)
{
    int m = 4, n = 5, k = 3;
    double A[k*m], B[k*n], C[m*n], D[m*n];

    A[0] = 1.0; A[3] = 2.0;  A[6] = 3.0; A[ 9] = 4.0;
    A[1] = 8.0; A[4] = 7.0;  A[7] = 6.0; A[10] = 5.0;
    A[2]  = 9.0;A[5] = 10.0; A[8] = 11.0;A[11] = 12.0;

    B[0] = 13.0; B[3] = 18.0; B[6] = 19.0; B[ 9] = 24.0; B[12] = 25.0;
    B[1] = 14.0; B[4] = 17.0; B[7] = 20.0; B[10] = 23.0; B[13] = 26.0;
    B[2] = 15.0; B[5] = 16.0; B[8] = 21.0; B[11] = 22.0; B[14] = 27.0;

    D[0] = 260.0; D[4] = 298.0; D[ 8] = 368.0; D[12] = 406.0; D[16] = 476.0;
    D[1] = 274.0; D[5] = 315.0; D[ 9] = 388.0; D[13] = 429.0; D[17] = 502.0;
    D[2] = 288.0; D[6] = 332.0; D[10] = 408.0; D[14] = 452.0; D[18] = 528.0;
    D[3] = 302.0; D[7] = 349.0; D[11] = 428.0; D[15] = 475.0; D[19] = 554.0;


    memset(C, 0, sizeof C);
    cblas_dgemm_scalar(CblasColMajor, CblasTrans, CblasNoTrans,
                       m, n, k, 1.0, A, k, B, k, 0.0, C, m);

    assert( memcmp(C, D, sizeof C) == 0 );
}

int main(int argc, char **argv)
{
    (void) argv;

    srand(time(NULL) + (long)&argc);
    test_matrix_print();
    test_matrix_allocate();
    test_matrix_one();
    test_vector_ddot();
    test_dgemm_scalar();

    bench_vector_ddot_incONE();
    bench_dgemm_scalar();

    return EXIT_SUCCESS;
}
