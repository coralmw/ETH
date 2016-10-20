#define N 2


#include <stdio.h>
#include <complex.h>
#include "kronecker_product.c"
//include "mkl.h"

const double Sz[N][N] = {[0][0] = 1, [0][1] = 0, [1][0] = 0, [1][1] = -1};
const complex Sy[N][N] = {[0][0] = 0, [0][1] = -1*I, [1][0] = 1*I, [1][1] = 0};
const double Sx[N][N] = {[0][0] = 0, [0][1] = 1, [1][0] = 1, [1][1] = 0};

int main() {

  double Hz[N*N][N*N], Hx[N*N][N*N], Hxz[N*N][N*N];
  Kronecker_Product(&Hz[0][0], (double *)Sz, N, N, (double*)Sz, N, N); 
  Kronecker_Product(&Hx[0][0], (double *)Sx, N, N, (double*)Sx, N, N); 
  Kronecker_Product(&Hxz[0][0], (double *)Sx, N, N, (double*)Sz, N, N); 

  complex Hy[N*N][N*N];      
  Kronecker_Product_complex(&Hy[0][0], (complex *)Sy, N, N, (complex*)Sy, N, N); // just use the complex version all the time

  printf("Combined hamiltonion for Sz of 2 spins:\n");
    for (int i = 0; i < N*N; i++) {
        for (int j = 0; j < N*N; j++) {
            printf("%4.1f,", Hz[i][j]);
        }
    printf("\n");
    }

  printf("Combined hamiltonion for Sx of 2 spins:\n");
    for (int i = 0; i < N*N; i++) {
        for (int j = 0; j < N*N; j++) {
            printf("%4.1f,", Hx[i][j]);
        }
    printf("\n");
    }

  printf("Combined hamiltonion for spin 1 mesured in Sz and 2 measured in Sx:\n");
    for (int i = 0; i < N*N; i++) {
        for (int j = 0; j < N*N; j++) {
            printf("%4.1f,", Hxz[i][j]);
        }
    printf("\n");
    }

  printf("Combined hamiltonion for Sy of 2 spins:\n");
    for (int i = 0; i < N*N; i++) {
        for (int j = 0; j < N*N; j++) {
            printf("%4.1f+%4.1fi,", crealf(Hy[i][j]), cimagf(Hy[i][j]));
        }
    printf("\n");
    }


  // feast params
  MKL_INT* fpm[64];
  feastinit (fpm);

  double *E, *res, *X;  // feast output params
  double epsout, trace;
  int i, loop, info, M;
  

  zfeast_heev('F', // full matrix
              N*N, // number of eigenvecs
              N*N, // other dim of H
              fpm;
              -10,
              10,
              N*N, // 2x no of evals
              &epsout, 
              &loop,

              
}
