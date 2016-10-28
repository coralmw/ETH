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
