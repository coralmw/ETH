/* C source code is found in dgemm_example.c */

#define min(x,y) (((x) < (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
//#include "kronecker_product.c"

void Kronecker_Product(double *C, double *A, int nrows, int ncols, 
                                               double *B, int mrows, int mcols) {}

int main()
{
    int sn = 2;
    double *Sz;
    double *H;

    printf (" Allocating memory for matrices aligned on 64-byte boundary for better \n"
            " performance \n\n");
    Sz = (double *)mkl_malloc( sn*sn*sizeof( double ), 64 );
    H = (double *)mkl_malloc( (sn*sn)*(sn*sn)*sizeof( double ), 64 );

    if (Sz == NULL || H == NULL) {
      printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
      mkl_free(Sz);
      mkl_free(H);
      return 1;
    }

    printf (" Intializing matrix data \n\n");


    Sz[0] = 1;
    Sz[1] = 2;
    Sz[2] = 3;
    Sz[3] = 4;

    for (int i = 0; i < (sn*2)*(sn*2); i++) {
        H[i] = (double)(0);
    }

    //cblas_dger(CblasRowMajor,
    //           sn*2, sn*2, 1.0, // size of returned array
    //           Sz, 1, Sz, 1,    // input vectors and strides
    //           H, sn*2);        // output
    printf("Calc the product\n");
    Kronecker_Product(H, Sz, 2, 2, Sz, 2, 2);

    for (int i = 0; i < (sn*2)*(sn*2); i++) {
        printf("%f,", H[i]);
        if ((i+1)%(sn*2) == 0) printf("\n");
    }


    printf ("\n Deallocating memory \n\n");
    mkl_free(Sz);
    mkl_free(H);

    printf (" Example completed. \n\n");
    return 0;
}

