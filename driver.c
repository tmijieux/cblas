#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <malloc.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <string.h>

#include "cblas.h"
#include "cblas/cblas.h"
#include "perf/perf.h"

#include "util.h"
#include "ddot.h"
#include "dgemm.h"

#define NB_ITER 1000

void test_matrix_print(void)
{
  double M[4 * 6] ALIGNED(32);
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

void test_ddot(cblas_ddot_t ddot)
{
  double X[] ALIGNED(32) = { 1, 2, 3, 4, 5, 6 };
  assert( DEQUAL( ddot(6, X, 1, X, 1), 91.0, 0.01) );

  double *Y = memalign(32, sizeof*Y*100);
  assert( ((long)Y & 31)  == 0 );

  for (int i = 0; i < 100; ++i)
    Y[i] = i+1.0;
  assert( DEQUAL( ddot(100, Y, 1, Y, 1), 338350, 0.01) );
  free(Y);
}

void bench_ddot(cblas_ddot_t ddot)
{
  int m = 50;
  while ( m < 1000000 ) {
    double *v1, *v2;
    v1 = tdp_vector_new(m);
    v2 = tdp_vector_new(m);
    tdp_vector_rand(m, 42., DBL_MAX, v1);
    tdp_vector_rand(m, -37.0, 500.0, v2);

    tdp_cache_garbage();
    perf_t p1, p2;
    perf(&p1);
    for (int i = 0; i < NB_ITER; ++i)
      ddot(m, v1, 1, v2, 1);
    perf(&p2);

    perf_diff(&p1, &p2);
    printf("m = %6d | ", m);
    uint64_t nb_op = 2 * m * NB_ITER;
    printf("%10g Mflops | time(µs) = ", perf_mflops(&p2, nb_op));
    perf_printmicro(&p2);

    free(v1); free(v2);
    m *= 1.25;
  }
}

void bench_daxpy(cblas_daxpy_t daxpy)
{
  int m = 50;
  while ( m < 1000000 ) {
    double *v1, *v2;
    v1 = tdp_vector_new(m);
    v2 = tdp_vector_new(m);
    tdp_vector_rand(m, 42., DBL_MAX, v1);
    tdp_vector_rand(m, -37.0, 500.0, v2);

    tdp_cache_garbage();
    perf_t p1, p2;
    perf(&p1);
    for (int i = 0; i < NB_ITER; ++i)
      daxpy(m, 1.0, v1, 1, v2, 1);
    perf(&p2);

    perf_diff(&p1, &p2);
    printf("m = %6d | ", m);
    uint64_t nb_op = 2L * m * NB_ITER;
    printf("%10g Mflops | time(µs) = ", perf_mflops(&p2, nb_op));
    perf_printmicro(&p2);

    free(v1); free(v2);
    m *= 1.25;
  }
}

void bench_dgemm(cblas_dgemm_t dgemm)
{
  for (int m = 50; m <= 1000; m += 20) {
    double *M1, *M2, *M3;
    M1 = tdp_avx256_aligned_matrix_new(m, m);
    M2 = tdp_avx256_aligned_matrix_new(m, m);
    M3 = tdp_avx256_aligned_matrix_new(m, m);

    int ld = UPPER_LD(m, 32);
    tdp_matrix_rand(ld, m, M1, 42., DBL_MAX);
    tdp_matrix_rand(ld, m, M2, -37.0, 500.0);

    tdp_cache_garbage();
    perf_t p1, p2;
    perf(&p1);
    dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
	  m, m, m, 1.0, M1, ld, M2, ld, 0.0, M3, ld);
    perf(&p2);

    perf_diff(&p1, &p2);
    printf("%d ", m);
    unsigned long nb_op = 2UL * CUBE(m);
    printf("%g ", perf_mflops(&p2, nb_op));
    perf_printmicro(&p2);
    free(M1); free(M2); free(M3);
  }
}

void test_dgemm(cblas_dgemm_t dgemm)
{
  int m = 4, n = 5, k = 3;
  double A[k*m] ALIGNED(32);
  double B[k*n] ALIGNED(32);
  double C[m*n] ALIGNED(32);
  double D[m*n] ALIGNED(32);

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

  for (int i = 0; i < m*n; ++i)
    C[i] = (double) i;

  dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
	m, n, k, 1.0, A, k, B, k, 0.0, C, m);

  assert( memcmp(C, D, sizeof C) == 0 );
}

void test_daxpy(cblas_daxpy_t daxpy)
{
  int m = 20;
  double *A, *B, *C;

  A = tdp_vector_new(m);
  B = tdp_vector_new(m);
  C = tdp_vector_new(m);

  tdp_vector_one(m, 33.0, A);
  tdp_vector_zero(m, B);

  tdp_vector_one(m, 99.0, C);

  daxpy(m, 3.0, A, 1, B, 1);
  assert( memcmp(B, C, m*sizeof C[0]) == 0 );
}

void test_square_dgemm(cblas_dgemm_t dgemm)
{
  int m = 3, n = 3, k = 3;
  double A[k*m] ALIGNED(32);
  double B[k*n] ALIGNED(32);
  double C[m*n] ALIGNED(32);
  double D[m*n] ALIGNED(32);

  A[0] = 1.0; A[3] = 2.0;  A[6] = 3.0;
  A[1] = 8.0; A[4] = 7.0;  A[7] = 6.0;
  A[2]  = 9.0;A[5] = 10.0; A[8] = 11.0;

  B[0] = 13.0; B[3] = 18.0; B[6] = 19.0;
  B[1] = 14.0; B[4] = 17.0; B[7] = 20.0;
  B[2] = 15.0; B[5] = 16.0; B[8] = 21.0;

  D[0] = 260.0; D[3] = 298.0; D[ 6] = 368.0;
  D[1] = 274.0; D[4] = 315.0; D[ 7] = 388.0;
  D[2] = 288.0; D[5] = 332.0; D[ 8] = 408.0;

  for (int i = 0; i < m*n; ++i)
    C[i] = (double) i;

  dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
	m, n, k, 1.0, A, k, B, k, 0.0, C, m);
  tdp_matrix_print(n,m,C,n,stdout);
  assert( memcmp(C, D, sizeof C) == 0 );
}

void test_big_square_dgemm(cblas_dgemm_t dgemm)
{
  int m = 35, n = 35, k = 35;
  double *A = tdp_matrix_new(m, n);
  double *B = tdp_matrix_new(m, n);
  double *C = tdp_matrix_new(m, n);
  double *D = tdp_matrix_new(m, n);

  tdp_matrix_one(m, n, 2.0, A, m);
  tdp_matrix_fill(m, n, 42.0, B, m);
  tdp_matrix_fill(m, n, 84.0, D, m);

  for (int i = 0; i < m*n; ++i)
    C[i] = (double) i;

  dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
	m, n, k, 1.0, A, k, B, k, 0.0, C, m);
  assert( memcmp(C, D, m*n*sizeof C[0]) == 0 );
}

void bench_dgemv(cblas_dgemv_t dgemv)
{
  for (int m = 100; m <= 1000; m += 50) {
    double *A, *X, *Y;
    A = tdp_avx256_aligned_matrix_new(m, m);
    X = tdp_vector_new(m);
    Y = tdp_vector_new(m);

    int ld = UPPER_LD(m, 32);
    tdp_matrix_rand(ld, m, A, 42., DBL_MAX);
    tdp_vector_rand(m, -37.0, 500.0, X);
    tdp_vector_rand(m, -37.0, 500.0, Y);

    tdp_cache_garbage();
    perf_t p1, p2;
    perf(&p1);
    dgemv(CblasColMajor, CblasTrans, m, m, 1.0, A, ld, X, 1, 0.0, Y, 1);
    perf(&p2);

    perf_diff(&p1, &p2);
    printf("m = %6d | ", m);
    uint64_t nb_op = 2 * SQUARE(m);
    printf("%8g Mflops | time(µs) = ", perf_mflops(&p2, nb_op));
    perf_printmicro(&p2);

    free(A); free(X); free(Y);
  }
}

void test_dgemv(cblas_dgemv_t dgemv)
{
  int m = 10;
  int ld = UPPER_LD(m, 32);
  double *A, *X, *Y, *Z;

  A = tdp_avx256_aligned_matrix_new(m, m);
  tdp_matrix_one(m, m, 3.0, A, ld);

  X = tdp_vector_new(m);
  Y = tdp_vector_new(m);
  Z = tdp_vector_new(m);

  tdp_vector_one(m, 5.0, X);
  tdp_vector_one(m, 7.0, Y);

  dgemv(CblasColMajor, CblasTrans,
	m, m, 1.0, A, ld, X, 1, 0.0, Y, 1);

  tdp_vector_one(m, 15.0, Z);
  assert ( memcmp(Z, Y, m*sizeof*Y) == 0 );

  tdp_vector_one(m, 7.0, Y);
  tdp_vector_one(m, 29.0, Z);
  dgemv(CblasColMajor, CblasTrans,
	m, m, 1.0, A, ld, X, 1, 2.0, Y, 1);

  assert ( memcmp(Z, Y, m*sizeof*Y) == 0 );

  free(A);
  free(X); free(Y); free(Z);
}


#define TEST(type, symbol)                              \
  do {							\
    printf("Testing \"%s\" %s.\n", #type, #symbol);	\
    test_##type(symbol);				\
  }while(0)

#define BENCH(type, symbol)                     \
  do {						\
    /* printf("\nBenching %s:\n", #symbol);  */	\
    printf("%s\n", #symbol);			\
    bench_##type(symbol);			\
  }while(0)


static void tests(void)
{
  test_matrix_print();
  test_matrix_allocate();
  test_matrix_one();

  TEST(daxpy, daxpy_basic);
  TEST(daxpy, daxpy_basic1);
  TEST(daxpy, daxpy_avx256);

  TEST(dgemv, dgemv_basic);
  TEST(dgemv, dgemv_avx);

  TEST(ddot, ddot_basic_Thomas);
  TEST(ddot, ddot_basic_Fatima_Zahra);
  TEST(ddot, ddot_avx_256);

  TEST(dgemm, dgemm_scalar_Fatima_Zahra);
  TEST(dgemm, dgemm_scalar_Thomas);
  TEST(dgemm, dgemm_scalar2_Thomas);
  TEST(dgemm, dgemm_i);
  TEST(dgemm, dgemm_j);
  TEST(dgemm, dgemm_k);
  TEST(dgemm, dgemm_fast_OMP);
  //TEST(square_dgemm, dgemm_block);
  TEST(big_square_dgemm, dgemm_block);

#ifdef USE_MKL
  TEST(ddot, cblas_ddot);
  TEST(dgemm, cblas_dgemm);
#endif
}

static void benches(void)
{
  BENCH(dgemm, dgemm_block);
  BENCH(dgemm, dgemm_block);
  BENCH(dgemm, dgemm_fast_OMP);

  BENCH(daxpy, daxpy_basic);
  BENCH(daxpy, daxpy_basic1);
  BENCH(daxpy, daxpy_avx256);

  BENCH(dgemv, dgemv_basic);
  BENCH(dgemv, dgemv_avx);

  BENCH(dgemm, dgemm_OMP);
  BENCH(dgemm, dgemm_fast_sequential);


#ifdef USE_MKL
  BENCH(dgemm, cblas_dgemm);
  BENCH(ddot, cblas_ddot);
#endif

  BENCH(ddot, ddot_avx_256);

  BENCH(ddot, ddot_basic_Thomas);
  BENCH(ddot, ddot_basic_Fatima_Zahra);


  BENCH(dgemm, dgemm_scalar_Fatima_Zahra);
  BENCH(dgemm, dgemm_scalar_Thomas);
  BENCH(dgemm, dgemm_scalar2_Thomas);
  BENCH(dgemm, dgemm_i);
  BENCH(dgemm, dgemm_j);
  BENCH(dgemm, dgemm_k);
}

int main(int argc, char **argv)
{
  (void) argv;
  srand(time(NULL) + (long)&argc);
  //tests();
  // benches();
  /* BENCH(dgemm, dgemm_i); */
  /* BENCH(dgemm, dgemm_j); */
  /* BENCH(dgemm, dgemm_k); */
  BENCH(dgemm, dgemm_block);
    
  return EXIT_SUCCESS;
}
