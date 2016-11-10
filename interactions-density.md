##finding the decomposition

start with the defenition of a state as before

$$
\ket{\psi} = \sum_{i j}C_{i j}\ket{S_{i}}\ket{B_{j}}
$$

with energy eigenvalues $\{ E_{ij} \}$

but now define the combined states with interactions as

$$
\ket{\psi} = \sum_{n}C_{n}\ket{N_n}
$$

with energy eigenvalues $\{ E_{n} \}$

We now need to find the representation of a arbitary bath state in the interaction included basis in order to perform the partial trace.

equating the representations of $\ket{\psi}$, and taking a single bath state existing so $C_{i j} = 0$ unless $i=z$ where $z$ is the index of the desired bath state,

$$
\sum_{i}C_{i z}\ket{S_{i}}\ket{B_{z}} = \sum_{n}C_{n}\ket{N_n}
$$

we observe

$$
\ket{\psi}\left( \ket{S_{i}}\ket{B_{j}} \right) = \ket{\psi}\left( \ket{B_{j}}\ket{S_{i}} \right)
$$

As $\ket{\psi}$, as the seperated basis allows each element of the system to act only on it's corrosponding element. We use this to remove $\ket{B_{j}}$ from the bracket.

$$
\ket{B_{z}}\left( \sum_{i}C_{i z}\ket{S_{i}} \right) = \sum_{n}C_{n}\ket{N_n}
$$

left multiply by $\sum_{i'}C^{\star}_{i' z}\bra{S_{i'}}$,

$$
\left(\sum_{i'}C^{\star}_{i z}\bra{S_{i'}} \right) \ket{B_{z}}\left( \sum_{i}C_{i z}\ket{S_{i}} \right) = \left(\sum_{i'}C^{\star}_{i' z}\bra{S_{i'}} \right) \sum_{n}C_{n}\ket{N_n}
$$

$$
\ket{B_{z}} \left(\sum_{ii'}C_{i z}C^{\star}_{i z}\bra{S_{i'}}\ket{S_{i}} \right) = \left(\sum_{i'}C^{\star}_{i' z}\bra{S_{i'}} \right) \sum_{n}C_{n}\ket{N_n}
$$

$$
\ket{B_{z}} \left(\sum_{ii'}C_{i z}C^{\star}_{i z}\delta_{ii'} \right) = \left(\sum_{i'}C^{\star}_{i' z}\bra{S_{i'}} \right) \sum_{n}C_{n}\ket{N_n}
$$

$$
\ket{B_{z}} \left(\sum_{i}C_{i z}C^{\star}_{i z} \right) = \sum_{i'n}C^{\star}_{i' z}C_{n}\bra{S_{i'}}\ket{N_n}
$$

$$
\ket{B_{z}} = \frac{\sum_{i'n}C^{\star}_{i' z}C_{n}\bra{S_{i'}}\ket{N_n}}{\sum_{i}C_{i z}C^{\star}_{i z}}
$$

##Density matrix

Introduce $\rho = \ket{\psi} \bra{\psi}$

$$
\rho = (\sum_{q}C_{q} \ket{N_q})(\sum_{q'}C^{\star}_{q'} \bra{N_q'})
$$

Find the partial trace,

$$
\sum_{z}\bra{B_z}\rho\ket{B_z} = \sum_{z} \left( \frac{\sum_{in}C_{i z}C^{\star}_{n}\ket{S_{i}}\bra{N_n}} {\sum_{i}C^{\star}_{i z}C_{i z}} \right) \left(\sum_{q} C_{q} \ket{N_q} \right)\left(\sum_{q'}C^{\star}_{q'} \bra{N_q'} \right) \left( \frac{\sum_{i'n'}C^{\star}_{i' z}C_{n'}\bra{S_{i'}}\ket{N_n'}}{\sum_{i}C_{i z}C^{\star}_{i z}} \right)
$$

$$
\sum_{z}\bra{B_z}\rho\ket{B_z} = \frac{1}{ \left( \sum_{i}C_{i z}C^{\star}_{i z} \right)^2} \sum_{z} \sum_{in} \sum_{qq'} \sum_{i'n'} C_{i z}C^{\star}_{n}C_{q}C^{\star}_{q'}C^{\star}_{i' z}C_{n'} \ket{S_{i}} \bra{N_n}\ket{N_q} \bra{N_q'}\ket{N_n'} \bra{S_{i'}}
$$

use orthonormality, $q'=n', q=n$,

$$
\sum_{z}\bra{B_z}\rho\ket{B_z} = \frac{1}{ \left( \sum_{i}C_{i z}C^{\star}_{i z} \right)^2} \sum_{z} \sum_{in} \sum_{i'n'} C_{i z}C^{\star}_{n}C_{n}C^{\star}_{n'}C^{\star}_{i' z}C_{n'} \ket{S_{i}} \bra{S_{i'}}
$$


and find the density element $\hat{\rho}^{\text{reduced}}_{S_1,S_2}$, applying orthonormaity we obtain $i=1,i'=2$

$$
\hat{\rho}^{\text{reduced}}_{S_1,S_2} = \sum_{nn'} \norm{C_{n}}^2 \norm{C_{n'}}^2 \sum_{z} \frac{C_{1 z}C^{\star}_{2 z}}{ \left( \sum_{i} \norm{C_{i z}}^2 \right)^2}   
$$

$$
\hat{\rho}^{\text{reduced}}_{S_1,S_2} = \sum_{n} \norm{C_{n}}^2 \sum_{n'}\norm{C_{n'}}^2 \sum_{z} \frac{C_{1 z}C^{\star}_{2 z}}{ \left( \sum_{i} \norm{C_{i z}}^2 \right)^2}   
$$

$$
\hat{\rho}^{\text{reduced}}_{S_1,S_2} =\sum_{z} \frac{C_{1 z}C^{\star}_{2 z}}{ \left( \sum_{i} \norm{C_{i z}}^2 \right)^2}   
$$


Both sums over













$$
