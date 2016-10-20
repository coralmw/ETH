
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

The Hamiltonion for a state mesured in a inconsistant basis, for example $S_{z2}\otimes S_{y1}$, is not diagonal.

A good primer is http://electron6.phys.utk.edu/qm1/modules/m10/twospin.htm.
