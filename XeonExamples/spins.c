#include <complex.h>
#include <math.h>
#include <assert.h>

#include "mkl.h"
#include "mkl_types.h"

#include "kronecker_product.h"
#include "tools.h"

#define SPINS 12


const MKL_Complex8 Sx[2*2] = {
  {0,0}, {1,0},
  {1,0}, {0,0}
};

const MKL_Complex8 Sy[2*2] = {
  {0,0}, {0,-1},
  {0,1}, {0,0}
};

const MKL_Complex8 Sz[2*2] = {
  {1,0}, {0,0},
  {0,0}, {-1,0}
};

const MKL_Complex8 idty[2*2] = {
  {1,0}, {0,0},
  {0,0}, {1,0}
};

const MKL_Complex8 *Sn[3] = {Sx, Sy, Sz};

int states;

void SNN(MKL_Complex8 *result, int spin, int dim){
  // result is a states*states matrix that needs the kronecter product of
  // kron( [pauli_dim if n=spin else idty for spin in range(SPINS)])
  // as we can only do one at once, we need a loop with SPINS-1 iterations
  // we use result to hold intermediate results as we know it's big enough
  // we start from the final matrix and work back.
  MKL_Complex8 *tempstart = malloc( states*states*sizeof( MKL_Complex8 ));
  MKL_Complex8 *tempend = malloc( states*states*sizeof( MKL_Complex8 ));
  if (tempstart == NULL || tempend == NULL) {
    exit(-1);
  }

  for (int ti=0; ti < states*states; ti++ ) {
    tempstart[ti] = (MKL_Complex8){0, 0};
    tempend[ti] = (MKL_Complex8){0, 0};
  }

  // construct the end of the product chain. make sure to include the pauli matirx if needed
  if (spin == SPINS-1) Kronecker_Product_MKL_Complex8(tempstart, idty, 2, 2, Sn[dim], 2, 2);
  if (spin == SPINS-2) Kronecker_Product_MKL_Complex8(tempstart, Sn[dim], 2, 2, idty, 2, 2);
  if (spin < SPINS-2) Kronecker_Product_MKL_Complex8(tempstart, idty, 2, 2, idty, 2, 2);

  int tempsize = 4;

  for (int ki=0; ki<SPINS-2; ki++){ // allredy donw the first 2
    //printf("kron product for iteration %d, tempsize is %d out of %d\n", ki, tempsize, states);
    // we need to do it from the outside first, and allocate a matrix to hold the temp result.

    if (spin == ki){
      Kronecker_Product_MKL_Complex8(tempend, Sn[dim], 2, 2, tempstart, tempsize, tempsize);
    } else {
      Kronecker_Product_MKL_Complex8(tempend, idty, 2, 2, tempstart, tempsize, tempsize);
    }

    tempsize = tempsize*2; // expand the working matrix
    assert(tempsize <= states);
    // copy tempend to tempstart to use next iteration, and zero tempend
    for (int ti=0; ti < states*states; ti++ ) {
      tempstart[ti] = tempend[ti];
      //tempend[ti] = (MKL_Complex8){0, 0};
    }
  } // end for kin in SPINS-3

  // return the result
  for (int ti=0; ti < states*states; ti++ ) {
    result[ti] = tempstart[ti];
  }

  free(tempstart);
  free(tempend);
}





int main() {
  printf( "Start SPINS\n" );
  states = pow(2,SPINS); // global scope

  MKL_Complex8 *Snn[SPINS][3];
  // [SPINS][3]; // pointers to S(1,2)(x,y,z)
  printf( "allocating memory and calculating submatrixcies\n" );

  // allocate space for all the S(1,2,3...)(x,y,z) matricies
  int spin, dim;
  #pragma omp parallel for private(spin,dim)
  for (spin=0; spin<SPINS; spin++){
    for (dim=0; dim<3; dim++){
      Snn[spin][dim] = malloc( states*states*sizeof( MKL_Complex8 ));
      assert(Snn[spin][dim] != NULL);

      SNN(Snn[spin][dim], spin, dim);

      // if (SPINS == 2){ // don't know how to do more spins yet. nested kron
      //   if (spin==0) {
      //     Kronecker_Product_MKL_Complex8(Snn[spin][dim], Sn[dim], 2, 2, idty, 2, 2);
      //   } else if (spin == 1) {
      //     Kronecker_Product_MKL_Complex8(Snn[spin][dim], idty, 2, 2, Sn[dim], 2, 2);
      //   }
      // } else {
      //   SNN(Snn[spin][dim], spin, dim);
      // }

      //print_cmatrix("", states, states, Snn[spin][dim], states);

    }
  }

  printf( "calculating H\n" );

  MKL_Complex8 *H = (MKL_Complex8 *)malloc( states*states*sizeof( MKL_Complex8 ));

  // MKL_Complex8 *spin1 = (MKL_Complex8 *)malloc( states*states*sizeof( MKL_Complex8 ));
  // MKL_Complex8 *spin2 = (MKL_Complex8 *)malloc( states*states*sizeof( MKL_Complex8 ));

  assert(H != NULL);
  for (int i = 0; i<states*states; i++) H[i] = (MKL_Complex8) {0, 0};
  int startspin;
  dim = 0;
  for (startspin = 0; startspin < SPINS; startspin++){
    for (dim = 0; dim<3; dim++){
      // we matrix-mult S(dim)1 by S(dim)2 and add the result to H
      //printf( "next part of H H\n" );
      // SNN(spin1, startspin, dim);
      // SNN(spin2, (startspin+1)%SPINS, dim);

      cblas_cgemm(CblasRowMajor, // Layout
        CblasNoTrans, // take the transpose of a?
        CblasNoTrans, // take the transpose of B?
        states, // rows of A, C
        states, // cols of B, C
        states, // clos A, rows B
        &(MKL_Complex8){1, 0},   // prefactor of the multiplication
        Snn[startspin][dim], //A Array, size lda* m.
        states,      // Specifies the leading dimension of a as declared in the calling (sub)program.
        Snn[(startspin+1)%SPINS][dim], // B Array, size ldb by k. Before entry the leading n-by-k part of the array b must contain the matrix B.
        states,      // Specifies the leading dimension of b as declared in the calling (sub)program.
        &(MKL_Complex8){1, 0}, // prefactor of addition
        H, // C
        states //Specifies the leading dimension of c as declared in the calling (sub)program.
        );

      // cblas_cgemm(CblasRowMajor, // Layout
      //   CblasNoTrans, // take the transpose of a?
      //   CblasNoTrans, // take the transpose of B?
      //   states, // rows of A, C
      //   states, // cols of B, C
      //   states, // clos A, rows B
      //   &(MKL_Complex8){1, 0},   // prefactor of the multiplication
      //   spin1, //A Array, size lda* m.
      //   states,      // Specifies the leading dimension of a as declared in the calling (sub)program.
      //   spin2, // B Array, size ldb by k. Before entry the leading n-by-k part of the array b must contain the matrix B.
      //   states,      // Specifies the leading dimension of b as declared in the calling (sub)program.
      //   &(MKL_Complex8){1, 0}, // prefactor of addition
      //   H, // C
      //   states //Specifies the leading dimension of c as declared in the calling (sub)program.
      //   );

      //cblas_csscal(states*states, 5.0, H, 1);
      //printf( "finished this mm for H\n" );
    }
  }


  printf( " SPINS Example Program Results\n" );
  MKL_Complex8 trace = {0,0};
  for (int ti = 0; ti < states; ti++){
    trace.real += H[ti+states*ti].real;
    trace.imag += H[ti+states*ti].imag;
  }
  printf("trace of H is %f+i%f", trace.real, trace.imag);
  //print_cmatrix_noimag( "H", states, states, H, states);

  // calculate eigenstates etc


}
