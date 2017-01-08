#include <complex.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "mkl.h"
#include "mkl_types.h"

#include "kronecker_product.h"
#include "tools.h"

#define NBATH 6
#define NSYS  2
#define PRINT 1
#define TSTEPS 1000
#define DT 0.01

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

void test_orth(const MKL_Complex8 *eigenvectors) {
  MKL_Complex8 *orth = malloc(states*states*sizeof(MKL_Complex8));
  for (int i=0; i<states; i++) {
    for (int j=0; j<states; j++) {
      cdotu(&orth[i*states+j], &states, &eigenvectors[i], &states, &eigenvectors[j], &states);
    }
  }
  print_cmatrix("orthonormality of eigenvectors:", states, states, orth, states);
}

void pauliposn(int n, int dim, int spins, MKL_Complex8 *result);
void HfromConnectivityMatrix(int size, float* J, MKL_Complex8 *result);

int main(int argc, char **argv) {
  printf( "Start SPINS\n" );
  states = pow(2,SPINS); // global scope

  if (PRINT) printf("open results files");
  FILE *results, *eigenfile;
  char *filename = argc >= 2 ? argv[1] : "results";
  char eignfilename[50];
  strcpy(eignfilename, argv[1]);
  strcat(eignfilename, "_eigenvalues");
  results = fopen(filename, "w");
  eigenfile = fopen(eignfilename, "w");

  if (results == NULL || eigenfile == NULL) {
    printf("file open fail!\n");
    exit(1);
  }


    printf( "calculating H\n" );


  printf( "saving eigen\n" );
  for (int i=0; i<pow(2,NSYS); i++) fprintf(eigenfile, "%f ", eigen[i].real);

  printf( "saving rhoreduced\n" );
  for (int i=0; i<pow(2,NSYS); i++) {
    for (int j=0; j<pow(2,NSYS); j++) {
      fprintf(results, "%d %d ", i, j);
      for (int tstep=0; tstep<TSTEPS; tstep++){
        if (rhoreduced[i*systemstates+j][tstep].real > 1.1) printf("wtf at t %d\n", tstep);
        fprintf(results, "%f ", rhoreduced[i*systemstates+j][tstep].real);
      }
      fprintf(results, "\n");
    }
  }

  // calculate_rhoreduced(basis_weights, eigen_right, eigen, 0.0, rhoreduced);
  //
  // print_fmatrix("rho-reduced:", 2, 2, rhoreduced, 2);
  //
  // FILE *upup, *downdown;
  // upup = fopen("results/upup", "w");
  // downdown = fopen("results/downdown", "w");
  //
  // float t = 0;
  // for (int ti = 0; ti<TSTEPS; ti++) {
  //   calculate_rhoreduced(basis_weights, eigen_right, eigen, t, rhoreduced);
  //   if (rhoreduced[0] > 1 || rhoreduced[3] > 1) print_fmatrix("What went wrong", 2, 2, rhoreduced, 2);
  //   fprintf(upup, "%f ", rhoreduced[0]);
  //   fprintf(downdown, "%f ", rhoreduced[3]);
  //   t += DT;
  // }
  // fclose(upup);
  // fclose(downdown);

}

void pauliposn(int n, int dim, int spins, MKL_Complex8 *result){
  // result is a spins*spins matrix that needs the kronecter product of
  // kron( [pauli_dim if n=spin else idty for spin in range(SPINS)])
  // as we can only do one at once, we need a loop with spins-1 iterations
  // we use result to hold intermediate results as we know it's big enough
  // we start from the final matrix and work back.
  int size = pow(2, spins);
  MKL_Complex8 *tempstart = malloc( size*sizeof( MKL_Complex8 ));
  if (tempstart == NULL) {
    exit(-1);
  }

  for (int ti=0; ti < size; ti++ ) {
    tempstart[ti] = (MKL_Complex8){0, 0};
    result[ti] = (MKL_Complex8){0, 0};
  }

  // construct the end of the product chain. make sure to include the pauli matirx if needed
  if (n == spins-1) Kronecker_Product_MKL_Complex8(tempstart, idty, 2, 2, Sn[dim], 2, 2);
  if (n == spins-2) Kronecker_Product_MKL_Complex8(tempstart, Sn[dim], 2, 2, idty, 2, 2);
  if (n <  spins-2) Kronecker_Product_MKL_Complex8(tempstart, idty, 2, 2, idty, 2, 2);

  int tempsize = 4;

  for (int ki=0; ki<spins-2; ki++){ // allredy donw the first 2
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
    for (int ti=0; ti < size; ti++ ) {
      tempstart[ti] = result[ti];
    }
  } // end for kin in spins-2

  // return the result
  for (int ti=0; ti < size; ti++ ) {
    result[ti] = tempstart[ti];
  }

  free(tempstart);
}

void HfromConnectivityMatrix(int spins, float* J, MKL_Complex8 *result){

  int size = pow(2, spins);
  MKL_Complex8 *left, *right;
  left = malloc(size*sizeof(MKL_Complex8));
  right = malloc(size*sizeof(MKL_Complex8));

  for (int i=0; i<spins; i++){
    for (int j = 0; j<spins; j++){
      for (int dim = 1; dim <= 3; dim++) {
        if (J[i*spins*j] != 0) {
          pauliposn(i, dim, spins, left);
          pauliposn(j, dim, spins, right);

          cblas_cgemm(CblasRowMajor, // Layout
            CblasNoTrans, // take the transpose of a?
            CblasNoTrans, // take the transpose of B?
            spins, // rows of A, C
            spins, // cols of B, C
            spins, // clos A, rows B
            &(MKL_Complex8){1, 0},   // prefactor of the multiplication
            left, //A Array, size lda* m.
            spins,      // Specifies the leading dimension of a as declared in the calling (sub)program.
            right, // B Array, size ldb by k. Before entry the leading n-by-k part of the array b must contain the matrix B.
            spins,      // Specifies the leading dimension of b as declared in the calling (sub)program.
            &J[i*spins*j], // prefactor of addition
            result, // C
            spins //Specifies the leading dimension of c as declared in the calling (sub)program.
            );
        }
      }
    }
  }
}
