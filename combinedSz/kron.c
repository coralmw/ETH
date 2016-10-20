// ma, na - rows/cols in A
// ma, mb - rows/cols in B
// A, B - input
// 

#include "mkl.h"

void Kron(const int ma, const int na, const int mb, const int nb, const double a, const double b, const double *A, const double *B, double *C) {
    int i, j, I, J, M, N;
    double c;
    M = ma*mb;
    N = na*nb;
    for (j=0; j<na; j++) {
        for (i=0; i<na; i++) {
            I = i*mb;
            J = j*nb;
            c = b*A[cm(i,j,ma)];
            assign_matrix_block(I, J, mb, nb, M, N, a, c, B, C);
        }
    }
}
