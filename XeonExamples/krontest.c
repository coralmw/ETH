#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "kronecker_product.h"
#include "tools.h"

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

void pauliposn(int n, int dim, int spins, MKL_Complex8 * temp, MKL_Complex8 *result){

  int size = pow(2, spins);
  printf("spins %i size %i\n", spins, size);

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
    assert(tempsize <= size);
    // copy tempend to tempstart to use next iteration, and zero tempend
    for (int ti=0; ti < size; ti++ ) {
      temp[ti] = result[ti];
    }
  } // end for kin in spins-2

  // return the result
  for (int ti=0; ti < size; ti++ ) {
    result[ti] = temp[ti];
  }
}



int main() {

  int spins = 3;

  MKL_Complex8 *left = malloc(pow(2, spins)*sizeof(MKL_Complex8));
  MKL_Complex8 *temp = malloc(pow(2, spins)*sizeof(MKL_Complex8));

  pauliposn(1, 2, spins, temp, left);
  print_cmatrix_noimag("Sx kron Sx", spins, spins, left, spins);
  free(left);
  free(temp);

}
