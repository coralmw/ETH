#include <complex.h>
#include <math.h>
#include <assert.h>

#include "mkl.h"
#include "mkl_types.h"

#include "kronecker_product.h"
#include "tools.h"

#define SPINS 3
#define PRINT 1


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
  if (tempstart == NULL) {
    exit(-1);
  }

  for (int ti=0; ti < states*states; ti++ ) {
    tempstart[ti] = (MKL_Complex8){0, 0};
    result[ti] = (MKL_Complex8){0, 0};
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
      Kronecker_Product_MKL_Complex8(result, Sn[dim], 2, 2, tempstart, tempsize, tempsize);
    } else {
      Kronecker_Product_MKL_Complex8(result, idty, 2, 2, tempstart, tempsize, tempsize);
    }

    tempsize = tempsize*2; // expand the working matrix
    assert(tempsize <= states);
    // copy tempend to tempstart to use next iteration, and zero tempend
    for (int ti=0; ti < states*states; ti++ ) {
      tempstart[ti] = result[ti];
    }
  } // end for kin in SPINS-3

  // return the result
  for (int ti=0; ti < states*states; ti++ ) {
    result[ti] = tempstart[ti];
  }

  free(tempstart);
}

void SeperatedState(int system, int bath, float* state) {

  // can Kron product with single value and it is identity
  float* temp = malloc(states*sizeof(float));
  assert(temp != NULL);

  int spinUp;

  float next[2];
  int size = 1;
  state[0] = 1;
  int statesconsidered;

  for (int bathspin = 0; bathspin<SPINS-1; bathspin++){
    spinUp = bath & (1 << bathspin);
    next[0] = spinUp==0 ? 1.0 : 0.0;
    next[1] = spinUp==0 ? 0.0 : 1.0;
    statesconsidered = pow(2, bathspin);

    Kronecker_Product(temp, next, 1, 2, state, 1, statesconsidered);
    for (int i=0; i<statesconsidered; i++) state[i] = temp[i];
  }

  // do the system state
  next[0] = system==0 ? 1.0 : 0.0;
  next[1] = system==0 ? 0.0 : 1.0;
  Kronecker_Product(state, next, 1, 2, temp, 1, pow(2, SPINS));
  free(temp);
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
      if (PRINT) printf("state %d dim %d\n", spin, dim);
      if (PRINT) print_cmatrix("", states, states, Snn[spin][dim], states );

    }
  }

  printf( "calculating H\n" );

  // Assign a initial state
  MKL_Complex8 PSI[states];
  for (int i = 0; i<states; i++) PSI[i] = (MKL_Complex8){0,0};
  PSI[0] = (MKL_Complex8){1,0};

  // Generate all the states
  float state[states];
  for (int system = 0; system<2; system++ ){
    for (int bath = 0; bath<states-2; bath++) {
      SeperatedState(system, bath, state); // calculate the vector
      if (PRINT) printf("state with system:%d bath%d\n", system, bath);
      if (PRINT) {
        for (int i=0;i<states;i++) printf("%f, ",state[i]);
        printf("\n");
      }
    }
  }

  MKL_Complex8 *H = (MKL_Complex8 *)malloc( states*states*sizeof( MKL_Complex8 ));

  assert(H != NULL);
  for (int i = 0; i<states*states; i++) H[i] = (MKL_Complex8) {0, 0};
  int startspin;
  dim = 0;
  for (startspin = 0; startspin < SPINS; startspin++){
    for (dim = 0; dim<3; dim++){
      // we matrix-mult S(dim)1 by S(dim)2 and add the result to H
      //printf( "next part of H H\n" );

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
      //printf( "finished this mm for H\n" );
    }
  }


  printf( " SPINS Example Program Results\n" );
  MKL_Complex8 trace = {0,0};
  for (int ti = 0; ti < states; ti++){
    trace.real += H[ti+states*ti].real;
    trace.imag += H[ti+states*ti].imag;
  }
  printf("trace of H is %f+i%f\n", trace.real, trace.imag);
  if (PRINT) print_cmatrix_noimag( "H", states, states, H, states);

  for (spin=0; spin<SPINS; spin++)
    for (dim=0; dim<3; dim++) free(Snn[spin][dim]); // remove the spins matrixes

  // find eigenvalues/vectors of H
  MKL_Complex8 *eigen = malloc( (states)*sizeof( MKL_Complex8 ));

  MKL_Complex8 *eigen_left = malloc( states*states*sizeof( MKL_Complex8 ));
  MKL_Complex8 *eigen_right = malloc( states*states*sizeof( MKL_Complex8 ));

  int info = LAPACKE_cgeev(
    LAPACK_ROW_MAJOR, //matrix Layout
    'N', // compute left eigenvectors
    'V',// compute right eigenvectors
    states, // size of H
    H, // matrix to eigensolve
    states, // lda
    eigen, // output for
    eigen_left, states, eigen_right, states
  );
   // H is nonsense after this call, but we don't need it
  /* Check for convergence */
  if( info > 0 ) {
    printf( "The algorithm failed to compute eigenvalues.\n" );
    exit( 1 );
  }
  for (int eval = 0; eval<states; eval++) printf("%d eigenvalue is %f+i%f\n", eval, eigen[eval].real, eigen[eval].imag);

  // create the eigenvalue matrix
  // MKL_Complex8 *eigenvectors = malloc( states*states*sizeof( MKL_Complex8 ));
  // for (int i=0; i<states; i++){
  //   for (int j = 0; j<states; j++) {
  //     eigenvectors[i+states*j].real = eigen_real[i+states*j];
  //     eigenvectors[i+states*j].imag = eigen_imag[i+states*j];
  //   }
  // }

  // for eigenvector in matrix, work out the dot of the PSI into it
  MKL_Complex8 *basis_weights = malloc( states*states*sizeof( MKL_Complex8 ));
  for (int i=0; i<states; i++){
    const int one = 1;
    cdotu(&basis_weights[i], &states, PSI, &one, &eigen_right[i], &states);
    printf("coeff of %d th basis state: %f\n", i, basis_weights[i].real);
  }


}
