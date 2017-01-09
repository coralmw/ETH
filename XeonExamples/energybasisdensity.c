#include <complex.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "mkl.h"
#include "mkl_types.h"

#include "kronecker_product.h"
#include "tools.h"

#define NBATH 2
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

void test_orth(int size, const float *eigenvectors) {
  float *orth = malloc(size*size*sizeof(float));
  for (int i=0; i<size; i++) {
    for (int j=0; j<size; j++) {
      orth[i*size+j] = sdot(&size, &eigenvectors[i], &size, &eigenvectors[j], &size);
    }
  }
  print_fmatrix("orthonormality of eigenvectors:", size, size, orth, size);
  free(orth);
}

void pauliposn(int n, int dim, int spins, MKL_Complex8 *result);
void HfromConnectivityMatrix(int size, float* J, float *result);
void eigensolveH(int size, float* H, float *Energy);
void reduceddensityEnergyBasis(float t, float *BathStates, float* SystemStates, float* SystemEnergies)

int main(int argc, char **argv) {
  printf( "Start SPINS\n" );
  states = pow(2,NSYS+NBATH); // global scope
  int matrixsize = pow(2, (NSYS+NBATH)*2);
  //
  // if (PRINT) printf("open results files");
  // FILE *results, *eigenfile;
  // char *filename = argc >= 2 ? argv[1] : "results";
  // char eignfilename[50];
  // strcpy(eignfilename, argv[1]);
  // strcat(eignfilename, "_eigenvalues");
  // results = fopen(filename, "w");
  // eigenfile = fopen(eignfilename, "w");
  //
  // if (results == NULL || eigenfile == NULL) {
  //   printf("file open fail!\n");
  //   exit(1);
  // }
  printf( "creating J\n" );

  float *Jsep, *Jcon;
  Jsep = malloc(pow(2, NSYS+NBATH) * sizeof(float));
  Jcon = malloc(pow(2, NSYS+NBATH) * sizeof(float)); // we also wuse this matrix for the BATH, SYS J
  if (Jsep == NULL || Jcon == NULL) exit(1);
  memset(Jsep, 0, sizeof Jsep);
  memset(Jcon, 0, sizeof Jcon);

  for (int i=0; i<NSYS+NBATH; i++) { // along the line above the diagonal
    if (i != NSYS-1) Jsep[i + (NSYS+NBATH)*(i+1)] = 1; // all but the system-bath link connected
  }

  for (int i=0; i<NSYS+NBATH-1; i++) { // along the line above the diagonal
    Jcon[i + (NSYS+NBATH)*(i+1)] = 1;  // all connected
  }
  Jcon[(NSYS+NBATH-1)*(NSYS+NBATH)] = 1; // connection from last bath to first system - topright.

  print_fmatrix("sep J matrix", NSYS+NBATH, NSYS+NBATH, Jsep, NSYS+NBATH);

  float *Hsep = malloc(matrixsize*sizeof(float));
  float *Hcon = malloc(matrixsize*sizeof(float));
  float *Hsys = malloc(pow(2, NSYS*2)*sizeof(float));
  float *Hbath = malloc(pow(2, NBATH*2)*sizeof(float));

  printf( "calculating H\n" );

  HfromConnectivityMatrix(NSYS+NBATH, Jsep, Hsep);
  HfromConnectivityMatrix(NSYS+NBATH, Jcon, Hcon);
  HfromConnectivityMatrix(NSYS, Jcon, Hsys); // The system and bath are  fully connected in the ring model
  HfromConnectivityMatrix(NBATH, Jcon, Hbath); // The system and bath are  fully connected in the ring model

  // print_fmatrix("Hsep", states, states, Hsep, states);

  printf( "calculating Eigensystems\n" );


  float *EigenEsep = malloc(matrixsize*sizeof(float));
  float *EigenEcon = malloc(matrixsize*sizeof(float));
  float *EigenEsys = malloc(pow(2, NSYS)*sizeof(float));
  float *EigenEbath = malloc(pow(2, NSYS)*sizeof(float));

  // eigenvectors overwrite the hamiltion passed in, but we don't need H anymore
  eigensolveH(NSYS+NBATH, Hsep, EigenEsep);
  eigensolveH(NSYS+NBATH, Hcon, EigenEcon);
  eigensolveH(NSYS, Hsys, EigenEsys);
  eigensolveH(NBATH, Hbath, EigenEbath);

  float *EigenVsep = Hsep;
  print_fmatrix("eigenvectors of Hsep", states, states, EigenVsep, states);
  test_orth(pow(2, NSYS+NBATH), EigenVsep);

  float *EigenVcon = Hcon;
  float *EigenVsys = Hsys;
  float *EigenVbath = Hbath;

  test_orth(pow(2, NSYS+NBATH), EigenVcon);

  MKL_Complex8 *Phi = malloc(pow(2, NSYS+NBATH)*sizeof(MKL_Complex8));
  memset(Phi, 0, sizeof(Phi));

  for (int i=0; i<pow(2, NSYS+NBATH); i++) Phi[i].real = EigenVsep[2 + i*states];
  print_cmatrix_noimag("Phi inital", pow(2, NSYS+NBATH), 1, Phi, pow(2, NSYS+NBATH));


}

void pauliposn(int n, int dim, int spins, MKL_Complex8 *result){
  // result is a spins*spins matrix that needs the kronecter product of
  // kron( [pauli_dim if n=spin else idty for spin in range(SPINS)])
  // as we can only do one at once, we need a loop with spins-1 iterations
  // we use result to hold intermediate results as we know it's big enough
  // we start from the final matrix and work back.

  int size = pow(2, spins*2);
  MKL_Complex8 *temp = malloc(size*sizeof(MKL_Complex8));

  // construct the end of the product chain. make sure to include the pauli matirx if needed
  if (n == spins-1) Kronecker_Product_MKL_Complex8(temp, idty, 2, 2, Sn[dim], 2, 2);
  if (n == spins-2) Kronecker_Product_MKL_Complex8(temp, Sn[dim], 2, 2, idty, 2, 2);
  if (n <  spins-2) Kronecker_Product_MKL_Complex8(temp, idty, 2, 2, idty, 2, 2);

  int tempsize = 4;

  for (int ki=0; ki<spins-2; ki++){ // allredy donw the first 2
    //printf("kron product for iteration %d, tempsize is %d out of %d\n", ki, tempsize, states);
    // we need to do it from the outside first, and allocate a matrix to hold the temp result.

    if (n == ki){
      Kronecker_Product_MKL_Complex8(result, Sn[dim], 2, 2, temp, tempsize, tempsize);
    } else {
      Kronecker_Product_MKL_Complex8(result, idty, 2, 2, temp, tempsize, tempsize);
    }

    tempsize = tempsize*2; // expand the working matrix
    assert(tempsize <= states);
    // copy tempend to tempstart to use next iteration, and zero tempend
    for (int ti=0; ti < size; ti++ ) {
      temp[ti] = result[ti];
    }
  } // end for kin in spins-2

  // return the result
  for (int ti=0; ti < size; ti++ ) {
    result[ti] = temp[ti];
  }
  free(temp);
}

void HfromConnectivityMatrix(int spins, float* J, float *result){

  int size = pow(2, spins*2);
  int dimsize = pow(2, spins);
  MKL_Complex8 *left, *right;

  left = malloc(size*sizeof(MKL_Complex8));
  right = malloc(size*sizeof(MKL_Complex8));

  // has to be complex here as the pauli matrix in the y-direction is complex
  MKL_Complex8 *HcomplexTemp = malloc(size*sizeof(MKL_Complex8));
  memset(HcomplexTemp, 0, sizeof(HcomplexTemp));

  if (left == NULL || right == NULL || HcomplexTemp == NULL) {
    printf("error alloc matrix\n");
    exit(2);
  }
  for (int dim = 0; dim < 3; dim++) {
    for (int i=0; i<spins; i++){
      memset(left, 0, sizeof(left));
      pauliposn(i, dim, spins, left);

      for (int j = 0; j<spins; j++){
          memset(right, 0, sizeof(right));
          pauliposn(j, dim, spins, right);

          if (J[i+spins*j] != 0) { // if cgemm gets a 0 prefactor for addition, it zeros your matrix.
            cblas_cgemm(CblasRowMajor, // Layout
              CblasNoTrans, // take the transpose of a?
              CblasNoTrans, // take the transpose of B?
              dimsize, // rows of A, C
              dimsize, // cols of B, C
              dimsize, // clos A, rows B
              &(MKL_Complex8){1, 0},   // prefactor of the multiplication
              left, //A Array, size lda* m.
              dimsize,      // Specifies the leading dimension of a as declared in the calling (sub)program.
              right, // B Array, size ldb by k. Before entry the leading n-by-k part of the array b must contain the matrix B.
              dimsize,      // Specifies the leading dimension of b as declared in the calling (sub)program.
              &(MKL_Complex8){J[i+spins*j], 0}, // prefactor of addition
              HcomplexTemp, // C
              dimsize //Specifies the leading dimension of c as declared in the calling (sub)program.
              );
          }
      }
    }
  }
  // cast to real as the hamiltionion is real eventually
  printf("copy\n");

  for (int i=0; i<size; i++) { result[i] = HcomplexTemp[i].real; }
  free(left);
  free(right);

  free(HcomplexTemp);
}

void eigensolveH(int size, float* H, float *Energy) {
  int states = pow(2, size);
  int info = LAPACKE_ssyevd( LAPACK_ROW_MAJOR, 'V', 'U', states, H, states, Energy );
  if( info > 0 ) {
    printf( "The algorithm failed to compute eigenvalues.\n" );
    exit( 1 );
  }
}
