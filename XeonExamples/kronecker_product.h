#include <complex.h>
#include "mkl_types.h"


void Kronecker_Product(double *C, double *A, int nrows, int ncols,
                                               double *B, int mrows, int mcols);

void Kronecker_Product_complex(complex *C, const complex *A, int nrows, int ncols,
                                              const complex *B, int mrows, int mcols);

void Kronecker_Product_MKL_Complex8(MKL_Complex8 *C, const MKL_Complex8 *A, int nrows, int ncols,
                                                   const MKL_Complex8 *B, int mrows, int mcols);
