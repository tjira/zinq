# Configuration–Interaction

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of the Configuration–Interaction (CI) method and how it is implemented in our framework. If you want to understand how we can describe electron correlation by mixing different electronic configurations together using a variational approach, this guide is written step-by-step for you.

---

## I. Wavefunction Representation

To understand the Configuration–Interaction method, we must first look at how the Hartree–Fock method represents a molecule. Hartree–Fock assumes that all electrons occupy a single configuration of molecular orbitals, described by a single Slater determinant. While this is a good starting point, it is an approximation that ignores the detailed ways in which electrons avoid one another (electron correlation). The Configuration–Interaction method improves on this by expressing the many-electron wavefunction as a linear combination of many different Slater determinants. These determinants represent different electronic configurations: the reference ground-state determinant, determinants where one electron has been excited to a virtual orbital, determinants where two electrons have been excited, and so on.

### 1. Variational Wavefunction Expansion

We write the Configuration–Interaction wavefunction expansion as

$$
| \Psi_{\text{CI}} \rangle = c_0 | \Phi_0 \rangle + \sum_i c_i | \Phi_i \rangle
$$

where $| \Phi_0 \rangle$ is the reference Hartree–Fock Slater determinant representing the mean-field ground state of the system, $| \Phi_i \rangle$ represents the excited Slater determinants where one or more electrons have been promoted from occupied molecular orbitals to virtual (unoccupied) molecular orbitals, and $c_0$ and $c_i$ are the variational coefficients that we want to determine. In a determinant-based Configuration–Interaction approach, we represent each Slater determinant in our code as an array of occupied spin-orbital indices. By comparing these arrays, we can determine the excitation level of each determinant relative to the ground state, which simplifies the evaluation of the interaction energies.

---

## II. Hamiltonian Evaluation and Diagonalization

### 2. Slater–Condon Rules

To find the coefficients that minimize the energy of our system, we must construct a large matrix representing the electronic Hamiltonian operator $\hat{H}$ in the basis of our Slater determinants. The elements of this matrix are denoted as $H_{ij}=\langle\Phi_i|\hat{H}|\Phi_j\rangle$. Because our molecular orbitals are orthonormal and the Hamiltonian operator only contains interactions involving at most two electrons, we can use the Slater–Condon rules to simplify these matrix elements. These rules reduce the complicated many-electron integrals to simple one- and two-electron molecular integrals based on how many spin-orbitals differ between the two determinants.

If the two determinants are identical, the matrix elements of a one-body operator $\hat{F}$ (which describes kinetic energy and nuclear attraction) and a two-body operator $\hat{G}$ (which describes electron-electron repulsion) are calculated as

$$
\langle\Phi_A|\hat{F}|\Phi_A\rangle=\sum_i\langle\phi_i|\hat{f}|\phi_i\rangle
$$

and

$$
\langle\Phi_A|\hat{G}|\Phi_A\rangle=\frac{1}{2}\sum_{i,j}\langle\phi_i\phi_j||\phi_i\phi_j\rangle
$$

where the sums run over all occupied spin-orbitals $\phi_i$ and $\phi_j$, and the double bar indicates an antisymmetrized two-electron integral combining the classical Coulomb repulsion and the quantum exchange interaction.

If the two determinants differ by exactly one spin-orbital (where determinant $A$ has orbital $\phi_p$ and determinant $B$ has orbital $\phi_r$ at that position), the matrix elements are calculated as

$$
\langle\Phi_A|\hat{F}|\Phi_B\rangle=\langle\phi_p|\hat{f}|\phi_r\rangle
$$

and

$$
\langle\Phi_A|\hat{G}|\Phi_B\rangle=\sum_i\langle\phi_p\phi_i||\phi_r\phi_i\rangle
$$

where the sum runs over the spin-orbitals common to both determinants, representing the transition interactions.

If the two determinants differ by exactly two spin-orbitals (where determinant $A$ has orbitals $\phi_p$ and $\phi_q$, and determinant $B$ has orbitals $\phi_r$ and $\phi_s$), the one-body matrix element is zero, and the two-body matrix element is calculated as

$$
\langle\Phi_A|\hat{G}|\Phi_B\rangle=\langle\phi_p\phi_q||\phi_r\phi_s\rangle
$$

which directly measures the interaction between the two excitations. If the determinants differ by three or more spin-orbitals, the matrix elements are exactly zero because the Hamiltonian operator only contains interactions between at most two electrons at a time. When evaluating these matrix elements in our code, we compare the occupied spin-orbital indices of the two determinants. Because the wavefunction must be antisymmetric, swapping the order of any two electrons changes the sign of the determinant. To account for this, the code counts the number of permutations needed to align the indices of the two determinants and applies a phase factor of $-1$ for each permutation.

### 3. Hamiltonian Diagonalization

Applying the variational principle to find the coefficients that minimize the energy leads to a matrix eigenvalue equation

$$
\mathbf{H}_{\text{CI}}\mathbf{C}_k=E_k\mathbf{C}_k
$$

where $\mathbf{H}_{\text{CI}}$ is the Hamiltonian matrix representation in the Slater determinant basis, $\mathbf{C}_k$ is the eigenvector containing the coefficients for state $k$, and $E_k$ is the corresponding electronic energy. The total energy of the state is the sum of the electronic energy $E_k$ and the classical nuclear repulsion energy. We solve this eigenvalue equation using symmetric matrix diagonalization methods from `linear_algebra.zig`.

---

## III. Nuclear Derivatives

### 4. Analytical Nuclear Gradient

To calculate the forces acting on the nuclei, we must find the derivative of the CI energy with respect to the nuclear coordinates. The analytical gradient of the energy of state $k$ with respect to a nuclear coordinate $x$ is given by

$$
\frac{dE_k}{dx}=\frac{dV_{\text{nuc}}}{dx}+\sum_{a,b}C_{ak}C_{bk}\frac{dH_{\text{CI},ab}}{dx}
$$

where $V_{\text{nuc}}$ is the nuclear repulsion energy, and $C_{ak}$ and $C_{bk}$ are the components of the CI eigenvector. Because the CI coefficients are variationally optimized by diagonalizing the Hamiltonian matrix, their derivatives with respect to the coordinates do not contribute to the energy derivative. However, the molecular orbitals themselves are optimized for the Hartree–Fock reference state, not the CI state, meaning they are not variational with respect to the CI energy. Therefore, we must explicitly calculate the derivatives of the molecular orbital coefficients by solving the Coupled-Perturbed Hartree–Fock (CPHF) equations to account for how the orbitals change when the nuclei move.
