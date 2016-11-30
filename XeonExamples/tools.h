#include <stdio.h>
#include <complex.h>
#include "mkl_types.h"


extern void print_matrix( char* desc, int m, int n, const complex* a, int lda );
void print_cmatrix( char* desc, int m, int n, const MKL_Complex8* a, int lda );

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, const complex* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %3.2f+i%3.2f", creal(a[i+j*lda]), cimag(a[i+j*lda]));
                printf( "\n" );
        }
}

void print_fmatrix( char* desc, int m, int n, const float* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %3.2f", creal(a[i+j*lda]));
                printf( "\n" );
        }
}

void print_cmatrix( char* desc, int m, int n, const MKL_Complex8* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %3.2f+i%3.2f", a[i+j*lda].real, a[i+j*lda].imag);
                printf( "\n" );
        }
}

void print_cmatrix_noimag( char* desc, int m, int n, const MKL_Complex8* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %3.2f", a[i+j*lda].real);
                printf( "\n" );
        }
}


MKL_Complex8 cmul(MKL_Complex8 a, MKL_Complex8 b) {
  return (MKL_Complex8){
                      (a.real*b.real) - (a.imag*b.imag),
                      (a.real*b.imag) + (a.imag*b.real)
                     };
}

float cmag(MKL_Complex8 a) {
  return pow(a.real, 2)+pow(a.imag, 2);
}


float cexp_MKL(MKL_Complex8 v) {
  complex vprime = v.real + I*v.imag;
  return cexp(vprime);
}
