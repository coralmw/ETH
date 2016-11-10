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

We are looking at a 2 spin systems, where the combined system is given by:

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

The Eigensystem of this is given by
$$
\left(
\begin{array}{cccc}
 -\frac{1}{4} \left(3 \hbar ^2\right) & \frac{\hbar ^2}{4} & \frac{\hbar ^2}{4} & \frac{\hbar ^2}{4} \\
 \{0,-1,1,0\} & \{0,0,0,1\} & \{0,1,1,0\} & \{1,0,0,0\} \\
\end{array}
\right)
$$

But this is only a eigensystem for S^2, not for S_z - for the simultainus mesurement of both magnitude and direction, the eigenvectors are the same (singlet and triplet) but the eigenvalues are 0 for the singlet state and $\hbar$ for the triplet states.

$$
S^{combined}_{z} = S_{1z} + S_{2z} = \hbar \left(
\begin{array}{cccc}
 1 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & -1 \\
\end{array}
\right)
$$

S_z has eigenvalues $-\hbar, 0, \hbar$ - (same as a spin-1 particle).

For a spin-up system coupled to a single spin-down bath, the probabilities of observing the system in a spin-up or spin-down state is given by the following $\rho^{reduced}$ matrix and plot.

$$
\rho^{reduced}_{\uparrow\downarrow} = \left(
\begin{array}{cc}
 \cos ^2\left(\frac{t}{2}\right) & -\frac{1}{2} i \sin (t) \\
 \frac{1}{2} i \sin (t) & \sin ^2\left(\frac{t}{2}\right) \\
\end{array}
\right)
$$
![probabilities of observing the system in a spin-up (blue) or spin-down (orange) state](combinedSz/densityMatrix.pdf){ width=50% }


For the state $\frac{1}{\sqrt{2}} \left(\ket{\uparrow}\ket{\downarrow}\right)_{S}\ket{\downarrow}_{B}$, the \rho^{reduced} was calculated with the following mathematica code.

~~~~~~~~~~~~~~Mathematica
Phi = 1/\[Sqrt]2 ({0, 1, 0, 0} + {0, 0, 0, 1})
PhiEigenBasis =
 Solve[{Phi == decomposition}, {\[Alpha], \[Beta], \[Gamma], \[Delta]}]
PhiEigenBasisoft = Psioft /. PhiEigenBasis
\[Rho] = Assuming[{t \[Element] Reals, \[HBar] \[Element] Reals},
  FullSimplify[
   ConjugateTranspose[PhiEigenBasisoft].PhiEigenBasisoft]]
\[Rho] // MatrixForm
PartialTrace =
 1/2 {0, 1, 0, 1}.\[Rho].{0, 1, 0, 1} +
  1/2 {1, 0, 1, 0}.\[Rho].{1, 0, 1, 0}
Subscript[\[Rho], reduced] =
 FullSimplify[
  Table[Table[
    i.\[Rho].j, {j, {{1, 1, 0, 0}, {0, 0, 1, 1}}}], {i, {{1, 1, 0,
      0}, {0, 0, 1, 1}}}]]
Subscript[\[Rho], reduced] // MatrixForm
~~~~~~~~~~~~~~

$$
\rho^{reduced}_{\frac{1}{\sqrt{2}} \left(\ket{\uparrow}\ket{\downarrow}\right)_{S}\ket{\downarrow}_{B}} = \left(
\begin{array}{cc}
 \frac{1}{4} (\cos (t)+1) & \frac{1}{4} (\cos (t)-2 i \sin (t)+1) \\
 \frac{1}{4} (\cos (t)+2 i \sin (t)+1) & \frac{1}{4} (5-3 \cos (t)) \\
\end{array}
\right)
$$

#Appendix
