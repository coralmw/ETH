Define the basis of system and bath

$$
\ket{\psi} = \sum_{i j}C_{i j}\ket{S_{i}}\ket{B_{j}}
$$

take the time evolution of the state

$$
\ket{\psi(t)} = \sum_{i j}C_{i j} e^{\frac{-i E_{i j} t}{\hbar}} \ket{S_{i}}\ket{B_{j}}
$$

Introduce $\rho(t) = \ket{\psi(t)} \bra{\psi(t)}$

$$
\rho(t) = (\sum_{i j}C_{i j} e^{\frac{-i E_{i j} t}{\hbar}} \ket{S_{i}}\ket{B_{j}})(\sum_{i' j'} C^{\star}_{i' j'} e^{\frac{i E_{i' j'} t}{\hbar}} \bra{B_{j'}}\bra{S_{i'}})
$$

Simplify by introducing the partial trace, $\sum_{n}\bra{B_n}\rho\ket{B_n}$, and using the presemed orthonormality of the basis to eliminate the sums over $j, j'$.

$$
\sum_{n}\bra{B_n}\rho\ket{B_n}
= \bra{B_n} (\sum_{i j}C_{i j} e^{\frac{-i E_{i j} t}/{\hbar}} \ket{S_{i}}\ket{B_{j}})(\sum_{i' j'}C^{\star}_{i' j'} e^{\frac{i E_{i' j'} t}/{\hbar}} \bra{B_{j'}}\bra{S_{i'}}) \ket{B_n}
$$
$$
= \sum_{n}\sum_{i i'} C_{i n} C^{\star}_{i' n} e^{\frac{-i (E_{i n}-E_{i' n})t}{\hbar}} \ket{S_{i}}\bra{S_{i'}}
$$

Now introduce the element of the reduced density matrix, $\hat{\rho}^{\text{reduced}}_{S_1,S_2}(t) = \ket{S_1}\sum_{n}\bra{B_n}\rho\ket{B_n}\bra{S_2}$

$$
\hat{\rho}^{\text{reduced}}_{S_1,S_2}(t) = \ket{S_1}\sum_{n}\bra{B_n}\rho\ket{B_n}\bra{S_2}
$$
$$
\hat{\rho}^{\text{reduced}}_{S_1,S_2}(t) = \bra{S_1}(\sum_{n}\sum_{i i'} C_{i n} C^{\star}_{i' n} e^{\frac{-i (E_{i n}-E_{i' n})t}{\hbar}} \ket{S_{i}}\bra{S_{i'}})\ket{S_2}
$$

use orthonormaity again to eliminate sums over $i,i'$

$$
\hat{\rho}^{\text{reduced}}_{S_1,S_2}(t) = \sum_{n} C_{1 n} C^{\star}_{2 n} e^{\frac{-i (E_{1 n}-E_{2 n})t}{\hbar}}
$$
