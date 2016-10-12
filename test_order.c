#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "cblas/cblas.h"

void cblas_dgemm_i(const enum CBLAS_ORDER Order,
		   const enum CBLAS_TRANSPOSE TransA,
		   const enum CBLAS_TRANSPOSE TransB,
		   const int M, const int N,const int K, 
		   const double alpha, 
		   const double *A,const int lda, 
		   const double *B, const int ldb,
		   const double beta, 
		   double *C, const int ldc){
  if(TransA==CblasNoTrans)  
      for (int i=0; i<M; i++)
	for (int j=0; j<N; j++) {
	  C[i+j*ldc]=0;
	  for (int k=0; k<K; k++)
	    C[i+j*ldc]=C[i+j*ldc]+A[i+k*lda]*B[k+j*ldb];
	}
    else if(TransA==CblasTrans)  
      for (int i=0; i<M; i++)
	for (int j=0; j<N; j++) {
	  C[i+j*ldc]=0;
	  for (int k=0; k<K; k++)
	    C[i+j*ldc]=C[i+j*ldc]+A[k+i*lda]*B[k+j*ldb];
	}
}


void cblas_dgemm_j(const enum CBLAS_ORDER Order,
		   const enum CBLAS_TRANSPOSE TransA,
		   const enum CBLAS_TRANSPOSE TransB,
		   const int M, const int N,const int K, 
		   const double alpha, 
		   const double *A,const int lda, 
		   const double *B, const int ldb,
		   const double beta, 
		   double *C, const int ldc){
  if(TransA==CblasNoTrans)  
      for (int j=0;j<N; j++)
	for (int i=0; i<M; i++) {
	  C[i+j*ldc]=0;
	  for (int k=0; k<K; k++)
	    C[i+j*ldc]=C[i+j*ldc]+A[i+k*lda]*B[k+j*ldb];
	}
    else if(TransA==CblasTrans)
      for (int j=0; j<N; j++)
	for (int i=0; i<M; i++){
	  C[i+j*ldc]=0;
	  for (int k=0; k<K; k++)
	    C[i+j*ldc]=C[i+j*ldc]+A[k+i*lda]*B[k+j*ldb];
	}
}
void cblas_dgemm_k(const enum CBLAS_ORDER Order,
		   const enum CBLAS_TRANSPOSE TransA,
		   const enum CBLAS_TRANSPOSE TransB,
		   const int M, const int N,const int K, 
		   const double alpha, 
		   const double *A,const int lda, 
		   const double *B, const int ldb,
		   const double beta, 
		   double *C, const int ldc){
  if(TransA==CblasNoTrans)
    for (int k=0; k<K; k++)
      for (int j=0;j<N; j++)
	for (int i=0; i<M; i++) {
	  C[i+j*ldc]=C[i+j*ldc]+A[i+k*lda]*B[k+j*ldb];
	}
    else if(TransA==CblasTrans)  
    for (int k=0; k<K; k++)
      for (int j=0;j<N; j++)
	for (int i=0; i<M; i++) {
	    C[i+j*ldc]=C[i+j*ldc]+A[k+i*lda]*B[k+j*ldb];
	}
}

