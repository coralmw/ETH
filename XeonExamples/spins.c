#include <complex.h>
#include <math.h>
#include <assert.h>

#include "mkl.h"
#include "mkl_types.h"

#include "kronecker_product.h"
#include "tools.h"

#define SPINS 6
#define NSYS  2
#define PRINT 1
#define TSTEPS 1000
#define DT 0.01

const int one = 1;

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


void calculate_eigensystem_fake2spin(MKL_Complex8 *eigen, MKL_Complex8 *eigen_left, MKL_Complex8 *eigen_right) {
  float paireigen[4] = {-6, 2, 2, 2};
  float pairvectors[4*4] = {
    0, 1, 0, 0,
    1/sqrt(2), 0, 1/sqrt(2), 0,
    -1/sqrt(2), 0, 1/sqrt(2), 0,
    0, 0, 0, 1
  };


  for (int i=0; i<4; i++) eigen[i].real = paireigen[i];
  for (int i=0; i<4*4; i++) eigen_right[i].real = pairvectors[i];
}


const MKL_Complex8 *Sn[3] = {Sx, Sy, Sz};

int states;

void SNN(MKL_Complex8 *result, int spin, int dim);

void SeperatedState(int system, int bath, MKL_Complex8* state) {
  int idx = 0;
  for (int s = 0; s < NSYS; s++) idx += pow(2, SPINS-s-1) * ((system & (1<<s)) != 0);

  for (int b = 0; b < SPINS-NSYS; b++) idx += pow(2, b) * ((bath & (1<<b)) != 0);
  for (int i=0; i<states; i++) state[i] = i == idx ? (MKL_Complex8){1.0,0.0} :  (MKL_Complex8){0.0,0.0};
}

void test_orth(const MKL_Complex8 *eigenvectors) {
  MKL_Complex8 *orth = malloc(states*states*sizeof(MKL_Complex8));
  for (int i=0; i<states; i++) {
    for (int j=0; j<states; j++) {
      cdotu(&orth[i*states+j], &states, &eigenvectors[i], &states, &eigenvectors[j], &states);
    }
  }
  print_cmatrix("orthonormality of eigenvectors:", states, states, orth, states);
}


void weights_matrix(const MKL_Complex8 *eigenvectors, const MKL_Complex8 *PSI, MKL_Complex8 *basis_weights);

void calculate_rhoreduced(const MKL_Complex8 *basis_weights, const MKL_Complex8 *eigenvectors, const MKL_Complex8 *eigenenergies, const float t, float *rhoreduced);
void calculate_rhoreduced_andsave(const MKL_Complex8 *basis_weights, const MKL_Complex8 *eigenvectors, const MKL_Complex8 *eigenenergies, MKL_Complex8 **rhoreduced);

void calculate_eigensystem(MKL_Complex8 *eigen, MKL_Complex8 *eigen_left, MKL_Complex8 *eigen_right);


int main() {
  printf( "Start SPINS\n" );
  states = pow(2,SPINS); // global scope

  printf( "calculating H\n" );

  // Assign a initial state
  MKL_Complex8 PSI[states];
  SeperatedState(1, 1, PSI);
  // MKL_Complex8 du[states]; SeperatedState(1,0,du);
  // MKL_Complex8 dd[states]; SeperatedState(1,1,dd);
  // MKL_Complex8 ddd[states]; SeperatedState(1,3,ddd);
  // MKL_Complex8 dddd[states]; SeperatedState(1,4,dddd);
  //
  //
  // // for (int i=0; i<states; i++) PSI[i] = add_MKL(du[i], cmul((MKL_Complex8){0,1}, dd[i]));
  // for (int i=0; i<states; i++) PSI[i] =
  //                                       add_MKL(du[i],
  //                                       add_MKL(dd[i],
  //                                       add_MKL(ddd[i],
  //                                               dddd[i]
  //                                       )));

  if (PRINT) print_cmatrix("PSI", states, 1, PSI, states);


  if (PRINT) {
    float norm = 0;
    for (int i = 0; i<states; i++) norm += cmag(PSI[i]);
    printf("norm^2 of state: %f\n", norm);
  }

  // find eigenvalues/vectors of H
  MKL_Complex8 *eigen = malloc( (states)*sizeof( MKL_Complex8 ));
  MKL_Complex8 *eigen_left = malloc( states*states*sizeof( MKL_Complex8 ));
  MKL_Complex8 *eigen_right = malloc( states*states*sizeof( MKL_Complex8 ));

  calculate_eigensystem(eigen, eigen_left, eigen_right);

  // test_orth(eigen_right);

  if (PRINT) {
    for (int eval = 0; eval<states; eval++) {
      printf("%d eigenvalue is %f+i%f ", eval, eigen[eval].real, eigen[eval].imag);
      printf("with eigenevctor: ");
      for (int i = 0; i<states; i++) printf("%3.2f+i%3.2f, ", eigen_right[i*states + eval].real, eigen_right[i*states + eval].imag);
      printf("\n");

      float norm = 0; // check the norm
      for (int i = 0; i<states; i++) norm += (eigen_right[i*states + eval].real * eigen_right[i*states + eval].real) + (eigen_right[i*states + eval].imag * eigen_right[i*states + eval].imag);
      printf("norm: %f\n", sqrt(norm));
    }
  }

  // for eigenvector in matrix, work out the dot of the PSI into it
  MKL_Complex8 *basis_weights = malloc( states*sizeof( MKL_Complex8 ));
  printf( "calculate weights\n" );
  weights_matrix(eigen_right, PSI, basis_weights);

  if (PRINT) {
    for (int i=0; i<states; i++) {
      printf("weight (%3.2f,%3.2f) for state: (", basis_weights[i].real, basis_weights[i].imag);
      for (int j=0; j<states; j++) printf(" %3.2f,", eigen_right[i + j*states].real);
      printf(")\n");
    }
  }
  //print_cmatrix("weights matrix:", states, 1, basis_weights, states);

  float basis_norm = 0;
  for (int i=0; i<states; i++) basis_norm += cmag(basis_weights[i]);
  if (PRINT) {
    printf("NORM of the basis (should be 1): %3.2f\n", sqrt(basis_norm));
  }

  // calculate the reduced density matrix at t=0.

  int systemstates = pow(2, NSYS);
  MKL_Complex8 *rhoreduced[systemstates*systemstates];
  for (int i=0; i<systemstates*systemstates; i++) {
    rhoreduced[i] = malloc(TSTEPS*sizeof(MKL_Complex8));
    for (int j=0; j<TSTEPS; j++) rhoreduced[i][j] = (MKL_Complex8){0,0};
    assert(rhoreduced[i] != NULL);
  }

  printf( "calculate rhoreduced for all t\n" );

  calculate_rhoreduced_andsave(basis_weights, eigen_right, eigen, rhoreduced);

  printf( "rhoreduced t=0\n" );
  for (int i=0; i<pow(2,NSYS); i++) {
    for (int j=0; j<pow(2,NSYS); j++) {
      printf("%3.2f, ", rhoreduced[i*systemstates+j][0].real);
    }
    printf("\n");
  }

  FILE *results;
  results = fopen("results", "w");

  if (results == NULL) {
    printf("file open fail!\n");
    exit(1);
  }

  printf( "saving rhoreduced\n" );


  for (int i=0; i<pow(2,NSYS); i++) {
    for (int j=0; j<pow(2,NSYS); j++) {
      fprintf(results, "%d %d ", i, j);
      for (int tstep=0; tstep<TSTEPS; tstep++){
        if (rhoreduced[i*systemstates+j][tstep].real > 1.1) printf("wtf at t %d", tstep);
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


void weights_matrix(const MKL_Complex8 *eigenvectors, const MKL_Complex8 *PSI, MKL_Complex8 *basis_weights) {
  for (int eval=0; eval<states; eval++){
    cdotc(&basis_weights[eval], &states, PSI, &one, &eigenvectors[eval], &states);
    if (PRINT>1 && basis_weights[eval].real > 0.0001) {
      printf("basis[%d] is %3.2f\n", eval,basis_weights[eval].real);
      printf("  eigen_right: "); for(int i=0;i<states;i++) printf("%3.2f, ", eigenvectors[eval + i*states].real); printf("\n");
    }
  }
}


void calculate_rhoreduced_andsave(const MKL_Complex8 *basis_weights, const MKL_Complex8 *eigenvectors, const MKL_Complex8 *eigenenergies, MKL_Complex8 **rhoreduced) {
  int numbath = pow(2,SPINS-NSYS);
  int numsystem = pow(2, NSYS);

  for (int i=0; i<numsystem; i++){ // elements of the rho-reduced matrix
    for (int j=0; j<numsystem; j++){
      // sum over the number of bath states to take the trace
      #pragma omp parallel for
      for (int n=0; n<numbath; n++) {
        MKL_Complex8 update, expfactor;
        MKL_Complex8 leftbraket, rightbraket; // dot products stored in here
        MKL_Complex8 leftstate[states], rightstate[states];
        SeperatedState(i, n, leftstate); // calculate the vector of the system in the ith state and the bath in the nth state
        SeperatedState(j, n, rightstate);

        // then sum over all the bath states
        for (int m=0; m<states; m++) {
          cdotc(&leftbraket,  &states,  leftstate,        &one,    &eigenvectors[m], &states);
          for (int mp=0; mp<states; mp++) {
            // cdotc takes the conj of the first vector
            cdotc(&rightbraket, &states,  &eigenvectors[mp], &states, rightstate,     &one);
            update = cmul(leftbraket, cmul(rightbraket, cmul(basis_weights[m], conj_MKL(basis_weights[mp]))));

            float t = 0;
            for (int tstep=0; tstep<TSTEPS; tstep++){
              expfactor = (MKL_Complex8){0, (eigenenergies[mp].real-eigenenergies[m].real)*t};
              #pragma omp atomic
              rhoreduced[i*numsystem+j][tstep].real += update.real * cexp_MKL(expfactor);
              #pragma omp atomic
              rhoreduced[i*numsystem+j][tstep].imag += update.imag * cexp_MKL(expfactor);
              t+=DT;
            }

            if (PRINT>2 && cmul(update, update).real > 0.00001 ) {
              printf("n:%d m:%d mp:%d\n", n, m, mp);
              printf("update to rhoreduced[%d][%d]+=%f\n", i,j, update.real);
            }
          }
        }
      }
    }
  }
}


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

void calculate_eigensystem(MKL_Complex8 *eigen, MKL_Complex8 *eigen_left, MKL_Complex8 *eigen_right) {

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
      //if (PRINT) printf("state %d dim %d\n", spin, dim);
      //if (PRINT) print_cmatrix("", states, states, Snn[spin][dim], states );

    }
  }

  MKL_Complex8 *H = malloc( states*states*sizeof( MKL_Complex8 ));
  // has to be complex as we construct it due to the SNN's being complex

  printf( "calculating H\n" );

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

  for (spin=0; spin<SPINS; spin++)
    for (dim=0; dim<3; dim++) free(Snn[spin][dim]); // remove the spins matrixes

  if (PRINT) {
    printf( " SPINS Example Program Results\n" );
    MKL_Complex8 trace = {0,0};
    for (int ti = 0; ti < states; ti++){
      trace.real += H[ti+states*ti].real;
      trace.imag += H[ti+states*ti].imag;
    }
    printf("trace of H is %f+i%f\n", trace.real, trace.imag);
     print_cmatrix_noimag( "H", states, states, H, states);
  }

  float *Hreal = malloc( states*states*sizeof( float ));
  float *eigenvalues = malloc( states*sizeof( float ));

  for (int i=0; i<states*states; i++) Hreal[i] = H[i].real; // cast to real
  free(H);

  printf( "eigensolve\n" );

  int info = LAPACKE_ssyevd( LAPACK_ROW_MAJOR, 'V', 'U', states, Hreal, states, eigenvalues );

  if( info > 0 ) {
    printf( "The algorithm failed to compute eigenvalues.\n" );
    exit( 1 );
  }

  for (int i=0; i<states*states; i++) eigen_right[i].real = Hreal[i]; // back to complex for return
  for (int i=0; i<states; i++) eigen[i].real = eigenvalues[i];

  free(Hreal);
  free(eigenvalues);
}
