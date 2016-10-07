#include "omp.h"
#include "mkl.h"
#include <stdio.h>

#define SIZE 1000

int main(int args, char *argv[]){

double *a, *b, *c;
a = new double [SIZE*SIZE];
b = new double [SIZE*SIZE];
c = new double [SIZE*SIZE];

double alpha=1, beta=1;
int m=SIZE, n=SIZE, k=SIZE, lda=SIZE, ldb=SIZE, ldc=SIZE, i=0, j=0;
char transa='n', transb='n';

for( i=0; i<SIZE; i++){
for( j=0; j<SIZE; j++){
a[i*SIZE+j]= (double)(i+j);
b[i*SIZE+j]= (double)(i*j);
c[i*SIZE+j]= (double)0;
}
}
cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

printf("row a c ");
for ( i=0;i<10;i++){
printf("%d: %f %f ", i, a[i*SIZE], c[i*SIZE]);
}

omp_set_num_threads(1);

for( i=0; i<SIZE; i++){
for( j=0; j<SIZE; j++){
a[i*SIZE+j]= (double)(i+j);
b[i*SIZE+j]= (double)(i*j);
c[i*SIZE+j]= (double)0;
}
}
cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

printf("row a c ");
for ( i=0;i<10;i++){
printf("%d: %f %f ", i, a[i*SIZE],
c[i*SIZE]);
}

omp_set_num_threads(2);
for( i=0; i<SIZE; i++){
for( j=0; j<SIZE; j++){
a[i*SIZE+j]= (double)(i+j);
b[i*SIZE+j]= (double)(i*j);
c[i*SIZE+j]= (double)0;
}
}
cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

printf("row a c ");
for ( i=0;i<10;i++){
printf("%d: %f %f ", i, a[i*SIZE],
c[i*SIZE]);
}

delete [] a;
delete [] b;
delete [] c;
}

