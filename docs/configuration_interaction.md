# Configuration–Interaction

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of the Configuration–Interaction (CI) method. If you want to understand how we can describe electron correlation by mixing different electronic configurations together using a variational approach, this guide is written step-by-step for you.

---

## I. Wavefunction Representation

To understand the Configuration–Interaction method, we must first look at how the Hartree–Fock method represents a molecule. Hartree–Fock assumes that all electrons occupy a single configuration of molecular orbitals, described by a single Slater determinant. While this is a good starting point, it is an approximation that ignores the detailed ways in which electrons avoid one another (electron correlation). The Configuration–Interaction method improves on this by expressing the many-electron wavefunction as a linear combination of many different Slater determinants. These determinants represent different electronic configurations: the reference ground-state determinant, determinants where one electron has been excited to a virtual orbital, determinants where two electrons have been excited, and so on.

### 1. Variational Wavefunction Expansion

We write the Configuration–Interaction wavefunction expansion for an electronic state $k$ as

$$
|\Psi_k\rangle=\sum_a c_{ak}|\Phi_a\rangle
$$

where $|\Phi_0\rangle$ is the reference Hartree–Fock Slater determinant representing the mean-field ground state of the system, $|\Phi_a\rangle$ represents the excited Slater determinants where one or more electrons have been promoted from occupied molecular orbitals to virtual molecular orbitals, and $c_{ak}$ are the variational coefficients that define state $k$. In a determinant-based Configuration–Interaction approach, each Slater determinant is specified by its occupied spin-orbital indices, from which the excitation level relative to the reference determinant is determined.

---

## II. Hamiltonian Evaluation and Diagonalization

### 2. Slater–Condon Rules

To find the coefficients that minimize the energy of the system, we construct the matrix representing the electronic Hamiltonian operator $\hat{H}$ in the basis of Slater determinants. The total electronic Hamiltonian is partitioned into one-electron and two-electron operators as

$$
\hat{H}=\sum_i\hat{h}^{\text{core}}(i)+\sum_{i<j}\hat{g}(i,j)
$$

where $\hat{h}^{\text{core}}$ describes the one-electron kinetic energy and nuclear attraction, and $\hat{g}(i,j)=\frac{1}{r_{ij}}$ represents the electron-electron Coulomb repulsion. The matrix elements in the determinant basis are denoted as $H_{\text{CI},ab}=\langle\Phi_a|\hat{H}|\Phi_b\rangle$. Because the molecular orbitals are orthonormal and the Hamiltonian contains at most two-body interactions, the Slater–Condon rules simplify these matrix elements to one- and two-electron molecular spin-orbital integrals based on the number of spin-orbitals differing between determinants $a$ and $b$.

If the two determinants are identical ($a=b$), the matrix element is the expectation value of the Hamiltonian evaluated as

$$
\langle\Phi_a|\hat{H}|\Phi_a\rangle=\sum_i\langle\phi_i|\hat{h}^{\text{core}}|\phi_i\rangle+\sum_{i<j}\langle\phi_i\phi_j||\phi_i\phi_j\rangle
$$

where the sums run over all occupied spin-orbitals $\phi_i$ and $\phi_j$, and $\langle\phi_i\phi_j||\phi_i\phi_j\rangle=\langle\phi_i\phi_j|\phi_i\phi_j\rangle-\langle\phi_i\phi_j|\phi_j\phi_i\rangle$ is the antisymmetrized two-electron integral combining classical Coulomb repulsion and quantum exchange.

If the two determinants differ by exactly one spin-orbital (where determinant $a$ contains spin-orbital $\phi_p$ and determinant $b$ contains spin-orbital $\phi_r$), the matrix element is evaluated as

$$
\langle\Phi_a|\hat{H}|\Phi_b\rangle=(-1)^{\sigma_{ab}}\left(\langle\phi_p|\hat{h}^{\text{core}}|\phi_r\rangle+\sum_i\langle\phi_p\phi_i||\phi_r\phi_i\rangle\right)
$$

where the sum runs over the spin-orbitals shared by both determinants, and $(-1)^{\sigma_{ab}}$ is a permutation phase factor determined by the number of orbital transpositions required to align the common spin-orbital occupations between the two determinants.

If the two determinants differ by exactly two spin-orbitals (where determinant $a$ contains spin-orbitals $\phi_p$ and $\phi_q$, and determinant $b$ contains spin-orbitals $\phi_r$ and $\phi_s$), the one-electron contribution vanishes, and the two-electron matrix element is evaluated as

$$
\langle\Phi_a|\hat{H}|\Phi_b\rangle=(-1)^{\sigma_{ab}}\langle\phi_p\phi_q||\phi_r\phi_s\rangle
$$

which directly couples the two excited configurations. If the determinants differ by three or more spin-orbitals, the Hamiltonian matrix element is identically zero because $\hat{H}$ contains only one- and two-body operators.

### 3. Hamiltonian Diagonalization

Applying the variational principle to optimize the expansion coefficients leads to the matrix eigenvalue equation

$$
\mathbf{H}_{\text{CI}}\mathbf{c}_k=E_k\mathbf{c}_k
$$

where $\mathbf{H}_{\text{CI}}$ is the Hamiltonian matrix representation in the Slater determinant basis, $\mathbf{c}_k$ is the eigenvector containing the configuration coefficients for state $k$, and $E_k$ is the corresponding electronic energy eigenvalue. The total energy of state $k$ is the sum of the electronic energy $E_k$ and the nuclear repulsion energy $V_{\text{nuc}}$.

---

## III. Nuclear Derivatives

### 4. Analytical Nuclear Gradient

To determine the forces acting on the nuclei, we calculate the derivative of the total energy of state $k$ with respect to the Cartesian nuclear coordinates $x$. Differentiating the variational energy expression yields the analytical nuclear gradient as

$$
\frac{dE_{k,\text{tot}}}{dx}=\frac{dV_{\text{nuc}}}{dx}+\sum_{a,b}c_{ak}c_{bk}\frac{dH_{\text{CI},ab}}{dx}
$$

where $V_{\text{nuc}}$ is the nuclear repulsion energy, and $c_{ak}$ and $c_{bk}$ are the components of the CI eigenvector $\mathbf{c}_k$. Because the CI coefficients are variationally optimized by diagonalizing the Hamiltonian matrix, their derivatives with respect to nuclear displacements vanish by the Hellmann–Feynman theorem. However, the molecular orbitals are optimized for the Hartree–Fock reference rather than the CI wavefunction, so the energy is not stationary with respect to changes in the molecular orbital coefficients, meaning their nuclear derivatives contribute explicitly to $\frac{dH_{\text{CI},ab}}{dx}$.

### 5. Differentiation via CPHF and Dual Numbers

Evaluating the derivative of each Hamiltonian matrix element $\frac{dH_{\text{CI},ab}}{dx}$ requires propagating nuclear derivatives through the atomic orbital integrals, the four-index molecular orbital transformation, and the Slater–Condon rules. To accomplish this, mean-field response theory is combined with forward-mode automatic differentiation using dual numbers.

First, the Coupled-Perturbed Hartree–Fock (CPHF) equations are solved to obtain the first-order response of the molecular orbitals to nuclear displacements. This yields the molecular orbital coefficient derivatives $\frac{d\mathbf{C}}{dx}$, which describe the geometric relaxation of the molecular orbital basis as the nuclei move.

Second, rather than manually differentiating the four-index molecular orbital transformation and the determinant evaluation rules, the derivative of the CI Hamiltonian matrix is evaluated using forward-mode automatic differentiation. Seeding the atomic orbital core Hamiltonian $\mathbf{H}^{\text{core}}$, the two-electron repulsion integrals $\mathbf{g}$, and the molecular orbital coefficients $\mathbf{C}$ with their exact first-order derivatives ($\frac{d\mathbf{H}^{\text{core}}}{dx}$, $\frac{d\mathbf{g}}{dx}$, and $\frac{d\mathbf{C}}{dx}$) automatically propagates the derivatives through the integral transformation and Slater–Condon evaluation using dual numbers. The resulting dual component of the Hamiltonian matrix directly provides $\frac{dH_{\text{CI},ab}}{dx}$, which is contracted with the CI eigenvectors to complete the analytical nuclear gradient.
