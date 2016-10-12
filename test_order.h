#ifndef TDP_PROD_MATRIX
#define TDP_PROD_MATRIX

#include "cblas/cblas.h"

void affiche(int m, int n, double* a,int lda);

void cblas_dgemm_scalar(const enum CBLAS_ORDER Order,
		   const enum CBLAS_TRANSPOSE TransA,
		   const enum CBLAS_TRANSPOSE TransB,
		   const int M, const int N,const int K, 
		   const double alpha, 
		   const double *A,const int lda, 
		   const double *B, const int ldb,
		   const double beta, 
		   double *C, const int ldc);


void cblas_dgemm_i(const enum CBLAS_ORDER Order,
		   const enum CBLAS_TRANSPOSE TransA,
		   const enum CBLAS_TRANSPOSE TransB,
		   const int M, const int N,const int K, 
		   const double alpha, 
		   const double *A,const int lda, 
		   const double *B, const int ldb,
		   const double beta, 
		   double *C, const int ldc);


void cblas_dgemm_j(const enum CBLAS_ORDER Order,
		   const enum CBLAS_TRANSPOSE TransA,
		   const enum CBLAS_TRANSPOSE TransB,
		   const int M, const int N,const int K, 
		   const double alpha, 
		   const double *A,const int lda, 
		   const double *B, const int ldb,
		   const double beta, 
		   double *C, const int ldc);

void cblas_dgemm_k(const enum CBLAS_ORDER Order,
		   const enum CBLAS_TRANSPOSE TransA,
		   const enum CBLAS_TRANSPOSE TransB,
		   const int M, const int N,const int K, 
		   const double alpha, 
		   const double *A,const int lda, 
		   const double *B, const int ldb,
		   const double beta, 
		   double *C, const int ldc);

#endif //TDP_PRODUCT_MATRIX
