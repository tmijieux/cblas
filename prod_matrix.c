#include<stdio.h>
#include<assert.h>
#include"cblas/cblas.h"


void cblas_dgemm_scalar(const enum CBLAS_ORDER Order,
			  const enum CBLAS_TRANSPOSE TransA,
			  const enum CBLAS_TRANSPOSE TransB,
			  const int M, const int N,const int K, 
			  const double alpha, 
			  const double *A,const int lda, 
			  const double *B, const int ldb,
			  const double beta, 
			  double *C, const int ldc){
  assert(Order==CblasColMajor);
  if(TransA==CblasNoTrans)  
      for (int i=0; i<M; i++)
	for (int j=0; j<N; j++) {    
	  for (int k=0; k<K; k++)
	    C[i+j*ldc]=C[i+j*ldc]+A[i+k*lda]*B[k+j*ldb];
	}
    else if(TransA==CblasTrans)  
      for (int i=0; i<M; i++)
	for (int j=0; j<N; j++) {    
	  for (int k=0; k<K; k++)
	    C[i+j*ldc]=C[i+j*ldc]+A[k+i*lda]*B[k+j*ldb];
	}
}


void cblas_dgemm(const enum CBLAS_ORDER Order,
			  const enum CBLAS_TRANSPOSE TransA,
			  const enum CBLAS_TRANSPOSE TransB,
			  const int M, const int N,const int K, 
			  const double alpha, 
			  const double *A,const int lda, 
			  const double *B, const int ldb,
			  const double beta, 
			  double *C, const int ldc){

  
  int m = M/3;
  double* A1, B1;
  for(int i=0; i<3;i++){
    /* tdp_submatrix(B,m,i*m,i*m*ldb,m,m,A1); */
    /* tdp_submatrix(A,m,i*m,i*m*lda,m,m,B1); */
    /*cblas_dgemm_scalar(Order,TransA,TransB,m,m,m,alpha,A1,lda,B1,ldb,beta,C,ldc); */
  }

}


/* double cblas_ddot(int N,const double *X, int incX,const double *Y,int incY){ */
/*   double Z=0; */
/*   for(int i=0;i<N;i++) */
/*     Z=Z+X[i*incX]*Y[i*incY]; */
/*   return Z; */

/* } */

void affiche(int m, int n, double* a,int lda){
   for(int i=0;i<m;i++){ 
     for(int j=0;j<n;j++){  
       printf("%g ",a[i+j*lda]); 
     } 
     printf("\n"); 
  }
}

