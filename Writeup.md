
---
title: Non-equilibrium dynamics and thermalisation in simple quantum systems
author:
-  Thomas Parks
-  Andrew Ho
abstract: |
  We aim to understand and model the behaviour of a simplified quantum system as it thermalises. We will investigate the behaviour of a spin-1/2 system using a computer simulation. The Eigenstate thermalization hypothesis (ETH) states that for a system prepared in some initial state where the expectation value of an observable $\hat{O}$ is far from that given by the microcanonical ensemble of this system, the expectation value of $\hat{O}$ will ultimately evolve in time to its value predicted by a microcanonical ensemble, without the invocation of any random processes. We shall simulate non-equlibrium quantum systems consisting of a several spin system coupled to a large bath and demonstrate this process.
bibliography: ETH.bib
link-citations: true
...


# The Eigenstate Thermalisation Hypothesis (ETH)

The time evolution of quantum systems is fully determined by the Schrodinger equation[@Schrodinger_1926]. Operationg on a known inital state, the future states of the system can be totally determined without any probabilistic quantaties energing, even if attempting to mesure this state will in general involve a probability distribuion over mesuremnet results.

This property of complete deteminisim stands in contrast to the conventional way to determine the time evolution of a complex classical system. In this formalistim, states are probability distributions over syatems that are at each instant locally thermally equlibulirated, and the future state of the system can only be determined by the forms of it's probability distribution.

The Eigenstate Thermalisation Hypothesis attempts to bridge this gap by showing how highly non-equlibirum quantum states can evolve to resemble equlibulirated states[@Srednicki_1994].

# The density matrix and limiting our observations of a system

In

# Derivation of the reduced density matrix for a system-bath state

Define the basis of system and bath in both the system and bath indexed states that are not a eigenbasis of the Hamiltonion and a singly indexed state of the eigenfunctions of H.

$$
\ket{\psi} = \sum_{i j}C_{i j}\ket{S_{i}B_{j}}
$$
$$
\ket{\psi} = \sum_{n}C_{n}\ket{n}
$$

Right multiply by $\bra{S_{i}B_{j}}$, and obtain

$$
C_{i j} = \sum_{n} C_{n}\bra{S_{i}B_{j}}\ket{n}
$$

Using the result for a element of the reduced density matrix,

$$
\hat{\rho}^{\text{reduced}}_{S_1,S_2} = \sum_{n} C_{1n}C^{\star}_{2n}
$$

And substuting the result for the seperated indexed coeffiecents in terms of the basis coefficents,

$$
\hat{\rho}^{\text{reduced}}_{S_1,S_2} = \sum_{n} \sum_{mm'} C_{m}\bra{S_{1}B_{n}}\ket{m}  C^{\star}_{m'}\bra{m'}\ket{S_{2}B_{n}}
$$

And we can apply the time evolution operator to this element as we know the eigenenergies in the eigenfunction basis.

$$
\hat{\rho}^{\text{reduced}}_{S_1,S_2} = \sum_{n} \sum_{mm'} \bra{S_{1}B_{n}}\ket{m}  \bra{m'}\ket{S_{2}B_{n}} C_{m}C^{\star}_{m'} e^{\frac{-i (E_{m}-E_{m'})t}{\hbar}}
$$






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

The Eigensystem of this is given by
$$
\left(
\begin{array}{cccc}
 -\frac{1}{4} \left(3 \hbar ^2\right) & \frac{\hbar ^2}{4} & \frac{\hbar ^2}{4} & \frac{\hbar ^2}{4} \\
 \{0,-1,1,0\} & \{0,0,0,1\} & \{0,1,1,0\} & \{1,0,0,0\} \\
\end{array}
\right)
$$


# Calculated results

## $\ket{\psi} = \ket{\uparrow \downarrow}$

$$
\hat{\rho}^{\text{reduced}}(t) = \left(
\begin{array}{cc}
 \sin ^2\left(\frac{t \hbar }{2}\right) & 0 \\
 0 & \cos ^2\left(\frac{t \hbar }{2}\right) \\
\end{array}
\right)
$$

And this state oscillates bwteen the up and down energy expectation values.

![Time evolution of the system portion of the $\ket{\uparrow \downarrow}$ state](combinedSz/graphics/upup.pdf){ width=70% }

## $\ket{\psi} = \frac{1}{\sqrt{2}} ( \ket{\uparrow \downarrow} +\ket{\downarrow \downarrow} )$

$$
\hat{\rho}^{\text{reduced}}(t) = \left(
\begin{array}{cc}
\frac{1}{4} (3-\cos (t \hbar )) & \frac{1}{4} (\cos (t \hbar )+i \sin (t \hbar )+1) \\
\frac{1}{4} (\cos (t \hbar )-i \sin (t \hbar )+1) & \frac{1}{4} (\cos (t \hbar )+1) \\
\end{array}
\right)
$$

And this state oscillates bwteen the up and down energy expectation values.

![Time evolution of the system portion of the $\frac{1}{\sqrt{2}} (\ket{\uparrow \downarrow} +\ket{\downarrow \downarrow})$ state](combinedSz/graphics/super-upd-downd.pdf){ width=70% }



# Performing calculations in Mathematica

By extending the definiation of the Hamiltonion to




~~~~~~~~Mathematica
NoSpins = 2;
NoStates = 2^NoSpins;
NoBasis = 2^(NoSpins - 1) - 1;
PauliPosN[n_, l_] :=
  Table[If[currn == n, \[HBar]/2 PauliMatrix[dim],
    IdentityMatrix[2]], {currn, l}];
Sdimn[dim_, NSpins_] :=
  Apply[Dot,
   Table[Nest[ArrayFlatten,
     Apply[Outer, {Times}~Join~PauliPosN[pos, NSpins]], NSpins], {pos,
      NSpins}]];

H = ArrayFlatten[Sum[Sdimn[dim, NoSpins],
     {dim, {1, 2, 3}}]];
SeperatedStates[s_, b_] := Apply[
   Join, Apply[
    Outer, {Times}~Join~
     Table[{BitGet[b, n], Boole[BitGet[b, n] == 0]}, {n, 0,
       NoSpins - 2}]~Join~{{Boole[s == 1], Boole[s == 0]}}], {0,
    NoSpins - 2}];

CombinedStates = Table[Normalize[evec], {evec, Eigenvectors[H]}];
CombinedE = Eigenvalues[H];
phi = 1/\[Sqrt]2 (SeperatedStates[1, 0] + SeperatedStates[0, 0]);
CombinedCoeffs = Table[Dot[phi, CombinedStates[[k]]], {k, NoStates}];
\[Rho][s1_, s2_] := ExpToTrig[\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(n = 0\), \(NoBasis\)]\ \(
\*UnderoverscriptBox[\(\[Sum]\), \(m = 1\), \(Length[
      CombinedStates]\)]\(
\*UnderoverscriptBox[\(\[Sum]\), \(p = 1\), \(Length[
       CombinedStates]\)]Dot[Conjugate[SeperatedStates[s1, n]],
       CombinedStates[\([\)\(m\)\(]\)]] Dot[
       Conjugate[CombinedStates[\([\)\(p\)\(]\)]],
       SeperatedStates[s2,
        n]] CombinedCoeffs[\([\)\(m\)\(]\)] Conjugate[
       CombinedCoeffs[\([\)\(p\)\(]\)]]
\*SuperscriptBox[\(E\), \(\(-I\) \((CombinedE[\([\)\(m\)\(]\)] -
          CombinedE[\([\)\(p\)\(]\)])\) t/\[HBar]\)]\)\)\)
  ]
p = Simplify[
   ArrayFlatten[Table[Table[\[Rho][s1, s2], {s1, 0, 1}], {s2, 0, 1}]]];
~~~~~~~~

# Introduction to the Xeon Phi

# References
