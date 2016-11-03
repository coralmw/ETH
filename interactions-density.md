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

left multiply by $\sum_{i'}C^{\star}_{i z}\bra{S_{i'}}$,

$$
\left(\sum_{i'}C^{\star}_{i z}\bra{S_{i'}} \right) \ket{B_{z}}\left( \sum_{i}C_{i z}\ket{S_{i}} \right) = \left(\sum_{i'}C^{\star}_{i z}\bra{S_{i'}} \right) \sum_{n}C_{n}\ket{N_n}
$$

$$
\ket{B_{z}} \left(\sum_{ii'}C_{i z}C^{\star}_{i z}\bra{S_{i'}}\ket{S_{i}} \right) = \left(\sum_{i'}C^{\star}_{i z}\bra{S_{i'}} \right) \sum_{n}C_{n}\ket{N_n}
$$
