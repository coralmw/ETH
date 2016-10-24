
#Project Proposal
##Non-equilibrium dynamics and thermalisation in simple quantum systems

To understand and model the behaviour of a simplified quantum system as it thermalises. We will investigate the behaviour of a spin-1/2 system using a coumpter simulation. The Eigenstate thermalization hypothesis (ETH) suggests that for a system prepared in some initial state where the expectation value of an observable $\hat{O}$ is far from that given by the microcanonical ensemble of this system, the expectation value of $\hat{O}$ will ultimately evolve in time to its value predicted by a microcanonical ensemble, without the invocation of any random processes. We shall simulate non-equlibrium quantum systems and hope to demonstrate this process.

#States

You can have a spin-1/2 system in the Sz basis, with up and down eigenstates going as

$$
S_{z} = \frac{\hbar}{2}\begin{pmatrix}
1 & 0\\
0 & -1
\end{pmatrix}
$$

With eigenstates

$$up=\binom{1}{0} , down=\binom{0}{1}$$

IF you add anouther spin state, you gain the combined systems:

$$
\ket{\uparrow\uparrow} = \begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix},
\ket{\uparrow\downarrow} = \begin{pmatrix} 0 \\ 1 \\ 0 \\ 0 \end{pmatrix},
\ket{\downarrow\uparrow} = \begin{pmatrix} 0 \\ 0 \\ 1 \\ 0 \end{pmatrix},
\ket{\downarrow\downarrow} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 1 \end{pmatrix}
$$

The combined $S^{tot}$ is given by the operation of $S_{z1}+S_{z2}$ only on the correct part of the state - so $S_{z1}$ sees the first spin, and $S_{z2}$ the second.

$$
S^{tot}\ket{\uparrow\uparrow} = (S_{z1}+S_{z2})(\ket{\uparrow}_{1}+\ket{\uparrow}_{2}) = S_{z1}\ket{\uparrow}_{1} + 0 + S_{z2}\ket{\uparrow}_{2} + 0 = \frac{\hbar}{2}\ket{\uparrow}_{1}+\frac{\hbar}{2}\ket{\uparrow}_{2} = \hbar\ket{\uparrow\uparrow}
$$

And the matrix form of the operator $S^{tot}$ is given by

$$
S^{tot} = \begin{bmatrix}
\bra{\uparrow\uparrow}S^{tot}\ket{\uparrow\uparrow} & \bra{\uparrow\uparrow}S^{tot}\ket{\uparrow\downarrow} & \bra{\uparrow\uparrow}S^{tot}\ket{\downarrow\uparrow} & \bra{\uparrow\uparrow}S^{tot}\ket{\downarrow\downarrow}\\
\bra{\uparrow\uparrow}S^{tot}\ket{\uparrow\uparrow} & \bra{\uparrow\downarrow}S^{tot}\ket{\uparrow\downarrow} & \bra{\downarrow\uparrow}S^{tot}\ket{\downarrow\uparrow} & \bra{\downarrow\downarrow}S^{tot}\ket{\downarrow\downarrow}\\
\bra{\uparrow\uparrow}S^{tot}\ket{\uparrow\uparrow} & \bra{\uparrow\downarrow}S^{tot}\ket{\uparrow\downarrow} & \bra{\downarrow\uparrow}S^{tot}\ket{\downarrow\uparrow} & \bra{\downarrow\downarrow}S^{tot}\ket{\downarrow\downarrow}\\
\bra{\uparrow\uparrow}S^{tot}\ket{\uparrow\uparrow} & \bra{\uparrow\downarrow}S^{tot}\ket{\uparrow\downarrow} & \bra{\downarrow\uparrow}S^{tot}\ket{\downarrow\uparrow} & \bra{\downarrow\downarrow}S^{tot}\ket{\downarrow\downarrow}\\
\end{bmatrix} = \frac{\hbar}{2}\begin{bmatrix}
2 & 0 & 0 & 0\\
0 & 0 & 0 & 0\\
0 & 0 & 0 & 0\\
0 & 0 & 0 & -2
\end{bmatrix}
$$

This matrix can be found by a summation over the Kronecker products of the basis matricies and the Identity matrix.

~~~~~~~~~Mathematica
ArrayFlatten[
  Outer[Times, PauliMatrix[3], IdentityMatrix[2]] +
   Outer[Times, IdentityMatrix[2], PauliMatrix[3]]
  ] // MatrixForm
~~~~~~~~~

Tha hamiliton of this system can be formed using the Kronecker product of the individual basis elements. This is implimented in Mathematica as Outer[Times, Sz, Sz], and a custom impl was used in C. This matrix is diagonal, and as such the eigenvalues are the diagonal elements and the eigenvectors are just unit vectors.

$$
H^{tot} = S_{z} \otimes S_{z} = \begin{bmatrix}
1 & 0 & 0 & 0\\
0 & -1 & 0 & 0\\
0 & 0 & -1 & 0\\
0 & 0 & 0 & 1
\end{bmatrix}
$$

The Eigensystem of this is given by
$$
\left(
\begin{array}{cccc}
 -\frac{1}{4} \left(3 \hbar ^2\right) & \frac{\hbar ^2}{4} & \frac{\hbar ^2}{4} & \frac{\hbar ^2}{4} \\
 \{0,-1,1,0\} & \{0,0,0,1\} & \{0,1,1,0\} & \{1,0,0,0\} \\
\end{array}
\right)
$$

The Hamiltonion for a state mesured in a inconsistant basis, for example $S_{z2}\otimes S_{y1}$, is not diagonal.

A good primer is http://electron6.phys.utk.edu/qm1/modules/m10/twospin.htm.

To find $S^2$, the total spin including interactions, you need to find $S_{1}\cdot S_{2}$. S is defined in the full 4 element basis.

$$
S_{1}\cdot S_{2} = S_{1_{x}} \cdot S_{2_{x}} + S_{1_{y}} \cdot S_{2_{y}} + S_{1_{z}} \cdot S_{2_{z}} = \frac{\hbar^2}{4} \left(
\begin{array}{cccc}
 1 & 0 & 0 & 0 \\
 0 & -1 & 2 & 0 \\
 0 & 2 & -1 & 0 \\
 0 & 0 & 0 & 1 \\
\end{array}
\right)
$$

#BLAS and LAPACK on the Xeon Phi

Intel provides the MKL - a coolection of numerical routienes that are optimised for use on intel processors. A kronecker product is not part of standerd BLAS/LAPACK so I used the following C code.

~~~~~~~~~~~~~C
void Kronecker_Product_complex(complex *C, complex *A, int nrows, int ncols,
                                               complex *B, int mrows, int mcols)
{
   int ccols, i, j, k, l;
   int block_increment;
   complex *pB;
   complex *pC, *p_C;

   ccols = ncols * mcols;
   block_increment = mrows * ccols;
   for (i = 0; i < nrows; C += block_increment, i++)
      for (p_C = C, j = 0; j < ncols; p_C += mcols, A++, j++)
         for (pC = p_C, pB = B, k = 0; k < mrows; pC += ccols, k++)
            for (l = 0; l < mcols; pB++, l++) *(pC+l) = *A * *pB;

}
~~~~~~~~~~~~~

#Meeting 2
#Density Matrix

$$
\hat{\rho}(t) = \ket{\psi(t)}\bra{\psi(t)} \tag*{Defn of density matrix}
$$
$$
\hat{\rho}^{\text{reduced}}_{i,j}(t) = \sum_{n} \ket{i}_{S}\ket{n}_{B} \hat{\rho} \bra{n}_{B}\bra{j}_{S}   \tag*{Defn of the reduced density matrix}
$$

The reduced density matrix, by the ETH, goes diagonal as the bath becomes suffucently complex.

$$
\ket{\psi} = \sum_{n} C_{n} \exp^{-\frac{i E_{n} t}{\hbar} } \ket{A_{n}}
$$

##Example 1

$$
\ket{\psi(0)} = ket{\uparrow}_{S}\ket{uparrow}_{B}
$$

##Solving $\hat{\rho}^{\text{reduced}}_{i,j}(t)$

$$
\hat{\rho}^{\text{reduced}}_{\uparrow,\uparrow}(t) = \ket{\uparrow}_{S}\ket{\uparrow}_{B}\hat{\rho}\bra{\uparrow}_{B}\bra{\uparrow}_{S} +
\ket{\uparrow}_{S}\ket{\downarrow}_{B}\hat{\rho}\bra{\downarrow}_{B}\bra{\uparrow}_{S}
$$

$$
\ket{\uparrow}_{B}\hat{\rho}\bra{\uparrow}_{B} = \ket{\uparrow}_{B}\ket{\psi}_{B}\bra{\psi}_{B}\bra{\uparrow}_{B} = \frac{1}{4}\left[ \ket{3}\bra{3} + e^{-i \hbar t} \ket{3}\bra{1} - e^{i \hbar t} \ket{1}\bra{1} + \ket{1}\bra{1} \right]
$$

$$
\ket{3}\bra{1} = \begin{pmatrix}0\\1\\1\\0\end{pmatrix}\begin{pmatrix}0&&-1&&1&&0\end{pmatrix} = \begin{pmatrix}0&&0&&0&&0\\0&&-1&&1&&0\\0&&-1&&1&&0\\0&&0&&0&&0\end{pmatrix}
$$

Outer product of matricies:

$$
\ket{1}\bra{1} = \begin{pmatrix}0\\-1\\1\\0\end{pmatrix}\begin{pmatrix}0&&-1&&1&&0\end{pmatrix} = \begin{pmatrix}0&&0&&0&&0\\0&&1&&-1&&0\\0&&-1&&1&&0\\0&&0&&0&&0\end{pmatrix}
$$
$$
\ket{3}\bra{3} = \begin{pmatrix}0\\1\\1\\0\end{pmatrix}\begin{pmatrix}0&&1&&1&&0\end{pmatrix} = \begin{pmatrix}0&&0&&0&&0\\0&&1&&1&&0\\0&&1&&1&&0\\0&&0&&0&&0\end{pmatrix}
$$
$$
\ket{1}\bra{3} = \begin{pmatrix}0\\-1\\1\\0\end{pmatrix}\begin{pmatrix}0&&1&&1&&0\end{pmatrix} = \begin{pmatrix}0&&0&&0&&0\\0&&-1&&-1&&0\\0&&1&&1&&0\\0&&0&&0&&0\end{pmatrix}
$$

Putting it back into the expression for $\ket{\uparrow}_{B}\hat{\rho}\bra{\uparrow}_{B}$

$$
\ket{\uparrow}_{B}\hat{\rho}\bra{\uparrow}_{B} = \begin{pmatrix}0&&1&&0&&0\end{pmatrix}\begin{pmatrix}0&&0&&0&&0\\
  0 && 2i \sin{\hbar t} && 2 \cos{\hbar t} && 0 \\
  0 && -2 \cos{\hbar t} && - 2i \sin{\hbar t} && 0 \\
0 && 0 && 0&& 0 \end{pmatrix}
\begin{pmatrix} 0\\1\\0\\0\end{pmatrix} = i \sin{\hbar t}
$$









\newpage
#Appendix
##Mathematica code for $S_{1}\cdot S_{2}$. {#sec:appendix_one}

~~~~~~~~~~~~~Mathematica
Subscript[S, 1 z] =
 ArrayFlatten[
  Outer[Times, \[HBar]/2 PauliMatrix[3],
   IdentityMatrix[2]]]; Subscript[S, 2 z] =
 ArrayFlatten[
  Outer[Times, IdentityMatrix[2], \[HBar]/2 PauliMatrix[3]]];
Subscript[S, 1 y] =
 ArrayFlatten[
  Outer[Times, \[HBar]/2 PauliMatrix[2],
   IdentityMatrix[2]]]; Subscript[S, 2 y] =
 ArrayFlatten[
  Outer[Times, IdentityMatrix[2], \[HBar]/2 PauliMatrix[2]]];
Subscript[S, 1 x] =
 ArrayFlatten[
  Outer[Times, \[HBar]/2 PauliMatrix[1],
   IdentityMatrix[2]]]; Subscript[S, 2 x] =
 ArrayFlatten[
  Outer[Times, IdentityMatrix[2], \[HBar]/2 PauliMatrix[1]]];
ArrayFlatten[Subscript[S, 1 z] + Subscript[S, 2 z]] // MatrixForm;
ArrayFlatten[Subscript[S, 1 y] + Subscript[S, 2 y]] // MatrixForm;
ArrayFlatten[Subscript[S, 1 z] + Subscript[S, 2 z]] // MatrixForm;
Subscript[S, 1 y].Subscript[S, 2 y] // MatrixForm
ArrayFlatten[
  ArrayFlatten[
   Subscript[S, 1 x].Subscript[S, 2 x] +
    Subscript[S, 1 y].Subscript[S, 2 y] +
    Subscript[S, 1 z].Subscript[S, 2 z]]] // MatrixForm
~~~~~~~~~~~~~

##First meeting results
### $H^{tot}$ in each basis

$$
H^{tot}_{x} = \frac{\hbar ^2}{4}\left(
\begin{array}{cccc}
0 & 0 & 0 & 1 \\
0 & 0 & 1 & 0 \\
0 & 1 & 0 & 0 \\
1 & 0 & 0 & 0 \\
\end{array}
\right),
H^{tot}_{y} = \frac{\hbar ^2}{4}\left(
\begin{array}{cccc}
0 & 0 & 0 & -1 \\
0 & 0 & 1 & 0 \\
0 & 1 & 0 & 0 \\
-1 & 0 & 0 & 0 \\
\end{array}
\right),
H^{tot}_{z} = \frac{\hbar ^2}{4}\left(
\begin{array}{cccc}
1 & 0 & 0 & 0 \\
0 & -1 & 0 & 0 \\
0 & 0 & -1 & 0 \\
0 & 0 & 0 & 1 \\
\end{array}
\right)
$$
