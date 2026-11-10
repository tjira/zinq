# Density Functional Theory

This document provides a clear, simple, and mathematically rigorous explanation of Density Functional Theory (DFT) within the Kohn–Sham framework. Building upon the foundational Hohenberg–Kohn existence and variational theorems, this guide explains how the interacting many-electron problem is mapped onto an auxiliary non-interacting reference system, how the Kohn–Sham potential and Roothaan–Kohn–Sham equations are derived, how exchange-correlation functionals are classified, and how spatial integration is evaluated on molecular grids.

---

## I. Hohenberg–Kohn Theorems and the Energy Functional

### 1. Hohenberg–Kohn Existence and Variational Theorems

In many-electron quantum mechanics, the electronic Hamiltonian operator $\hat{H}$ is composed of the electronic kinetic energy $\hat{T}$, the electron-electron Coulomb repulsion $\hat{V}_{\text{ee}}$, and the interaction of electrons with an external nuclear electrostatic potential $v(\mathbf{r})$, expressed in terms of the electron density operator $\hat{\rho}(\mathbf{r})$ as

$$
\hat{H}=\hat{T}+\hat{V}_{\text{ee}}+\int v(\mathbf{r})\hat{\rho}(\mathbf{r})d\mathbf{r}
$$

where $v(\mathbf{r})=-\sum_A\frac{Z_A}{|\mathbf{r}-\mathbf{R}_A|}$ is the external nuclear attraction potential defining the one-electron core Hamiltonian operator $\hat{h}=-\frac{1}{2}\nabla^2+v(\mathbf{r})$. The first Hohenberg–Kohn theorem proves that the ground-state electronic density $\rho(\mathbf{r})=\langle\Psi|\hat{\rho}(\mathbf{r})|\Psi\rangle$ uniquely determines the external potential $v(\mathbf{r})$ up to a trivial additive constant. To demonstrate this by contradiction, assume that two distinct external potentials $v_1(\mathbf{r})$ and $v_2(\mathbf{r})$ differing by more than a constant yield non-degenerate ground-state wavefunctions $|\Psi_1\rangle$ and $|\Psi_2\rangle$ that share the exact same ground-state density $\rho(\mathbf{r})$. Applying the Rayleigh–Ritz variational inequality to the ground state of $\hat{H}_1$ gives

$$
E_1<\langle\Psi_2|\hat{H}_1|\Psi_2\rangle=\langle\Psi_2|\hat{H}_2|\Psi_2\rangle+\int\left(v_1(\mathbf{r})-v_2(\mathbf{r})\right)\rho(\mathbf{r})d\mathbf{r}=E_2+\int\left(v_1(\mathbf{r})-v_2(\mathbf{r})\right)\rho(\mathbf{r})d\mathbf{r}
$$

Interchanging the indices and evaluating the inequality for the ground state of $\hat{H}_2$ yields

$$
E_2<\langle\Psi_1|\hat{H}_2|\Psi_1\rangle=E_1+\int\left(v_2(\mathbf{r})-v_1(\mathbf{r})\right)\rho(\mathbf{r})d\mathbf{r}
$$

Adding these two strict inequalities produces the manifest contradiction $E_1+E_2<E_1+E_2$. Consequently, two distinct potentials cannot generate the same ground-state electron density, proving that the ground-state density $\rho(\mathbf{r})$ uniquely specifies the external potential $v(\mathbf{r})$, the Hamiltonian $\hat{H}$, and all ground- and excited-state electronic properties of the system.

The second Hohenberg–Kohn theorem defines a universal functional $F_{\text{HK}}[\rho]=\langle\Psi[\rho]|(\hat{T}+\hat{V}_{\text{ee}})|\Psi[\rho]\rangle$ that depends only on the density and is completely independent of the external potential. The total ground-state energy functional is expressed as

$$
E_v[\rho]=F_{\text{HK}}[\rho]+\int v(\mathbf{r})\rho(\mathbf{r})d\mathbf{r}
$$

For any trial density $\tilde{\rho}(\mathbf{r})$ associated with an external potential and integrating to the correct number of electrons $N$, the energy functional satisfies the variational condition $E_0\le E_v[\tilde{\rho}]$, where the minimum is attained if and only if the trial density matches the exact ground-state density.

### 2. The Kohn–Sham Auxiliary System

Although the Hohenberg–Kohn theorems establish density as the fundamental physical variable, the kinetic energy functional $T[\rho]$ of an interacting many-electron system cannot be expressed accurately as an explicit local functional of the density. Kohn and Sham solved this limitation by mapping the physical interacting system onto an auxiliary reference system of $N$ non-interacting electrons moving in an effective local single-particle potential $v_{\text{eff}}(\mathbf{r})$.

The ground state of this non-interacting system is represented by a single Slater determinant of $N$ orthonormal single-particle orbitals $\{\phi_i(\mathbf{r})\}$, whose exact non-interacting kinetic energy is computed analytically as

$$
T_{\text{s}}[\rho]=-\frac{1}{2}\sum_i^{\text{occ}}\int\phi_i^*(\mathbf{r})\nabla^2\phi_i(\mathbf{r})d\mathbf{r}
$$

The single-particle orbitals are chosen such that their sum of squared moduli reproduces the exact ground-state electronic density of the interacting system as

$$
\rho(\mathbf{r})=\sum_i^{\text{occ}}|\phi_i(\mathbf{r})|^2
$$

The universal functional $F_{\text{HK}}[\rho]$ is then partitioned by separating the non-interacting kinetic energy $T_{\text{s}}[\rho]$ and the classical Hartree electrostatic repulsion $J[\rho]=\frac{1}{2}\int\int\frac{\rho(\mathbf{r})\rho(\mathbf{r}')}{|\mathbf{r}-\mathbf{r}'|}d\mathbf{r}d\mathbf{r}'$, leaving the exchange-correlation energy functional $E_{\text{xc}}[\rho]$ defined as

$$
E_{\text{xc}}[\rho]=\left(T[\rho]-T_{\text{s}}[\rho]\right)+\left(V_{\text{ee}}[\rho]-J[\rho]\right)
$$

which encapsulates the non-classical quantum exchange, dynamical electron correlation, and the difference between the true interacting kinetic energy and the non-interacting single-particle kinetic energy.

---

## II. Kohn–Sham Equations and Basis Set Expansion

### 3. Variational Derivation of the Kohn–Sham Potential

The total electronic energy functional in the Kohn–Sham formulation is written in terms of the occupied single-particle orbitals as

$$
E_{\text{KS}}[\{\phi_i\}]=-\frac{1}{2}\sum_i^{\text{occ}}\int\phi_i^*(\mathbf{r})\nabla^2\phi_i(\mathbf{r})d\mathbf{r}+\int v(\mathbf{r})\rho(\mathbf{r})d\mathbf{r}+\frac{1}{2}\int\int\frac{\rho(\mathbf{r})\rho(\mathbf{r}')}{|\mathbf{r}-\mathbf{r}'|}d\mathbf{r}d\mathbf{r}'+E_{\text{xc}}[\rho]
$$

We minimize this energy functional with respect to orbital variations $\phi_i^*(\mathbf{r})$ subject to the orbital orthonormality constraint $\int\phi_i^*(\mathbf{r})\phi_j(\mathbf{r})d\mathbf{r}=\delta_{ij}$ using Lagrange multipliers $\epsilon_{ij}$. Setting the functional derivative to zero and diagonalizing the Hermitian multiplier matrix yields the canonical Kohn–Sham eigenvalue equations

$$
\hat{h}_{\text{KS}}\phi_i(\mathbf{r})=\left(-\frac{1}{2}\nabla^2+v_{\text{eff}}(\mathbf{r})\right)\phi_i(\mathbf{r})=\epsilon_i\phi_i(\mathbf{r})
$$

where $\epsilon_i$ are the Kohn–Sham orbital energies, and the effective local potential $v_{\text{eff}}(\mathbf{r})$ is defined as

$$
v_{\text{eff}}(\mathbf{r})=v(\mathbf{r})+v_{\text{H}}(\mathbf{r})+v_{\text{xc}}(\mathbf{r})
$$

In this effective potential, $v(\mathbf{r})=-\sum_A\frac{Z_A}{|\mathbf{r}-\mathbf{R}_A|}$ is the nuclear attraction potential, $v_{\text{H}}(\mathbf{r})=\int\frac{\rho(\mathbf{r}')}{|\mathbf{r}-\mathbf{r}'|}d\mathbf{r}'$ is the classical Hartree electrostatic potential, and $v_{\text{xc}}(\mathbf{r})$ is the local exchange-correlation potential defined by the functional derivative

$$
v_{\text{xc}}(\mathbf{r})=\frac{\delta E_{\text{xc}}[\rho]}{\delta\rho(\mathbf{r})}
$$

### 4. Roothaan–Kohn–Sham Matrix Equations

To solve the Kohn–Sham equations computationally for molecular systems, each single-particle orbital is expanded in a known atomic orbital basis set as $|\phi_i\rangle=\sum_\mu C_{\mu i}|\chi_\mu\rangle$. Projecting the differential eigenvalue equations onto the atomic basis functions $\langle\chi_\mu|$ transforms the problem into the algebraic Roothaan–Kohn–Sham generalized eigenvalue equations

$$
\mathbf{F}\mathbf{C}=\mathbf{S}\mathbf{C}\mathbf{E}
$$

where $\mathbf{S}$ is the atomic orbital overlap matrix with elements $S_{\mu\nu}=\langle\chi_\mu|\chi_\nu\rangle$, $\mathbf{C}$ is the molecular orbital expansion coefficient matrix, and $\mathbf{E}=\text{diag}(\epsilon_1,\epsilon_2,\dots,\epsilon_{N_{\text{AO}}})$ is the diagonal matrix of orbital energies. The Kohn–Sham Fock matrix elements in the atomic basis are evaluated as

$$
F_{\mu\nu}=H_{\mu\nu}+J_{\mu\nu}-c_{\text{x}}a_{\text{x}}K_{\mu\nu}+V_{\mu\nu}^{\text{xc}}
$$

where $H_{\mu\nu}=\langle\chi_\mu|\hat{h}|\chi_\nu\rangle$ is the core Hamiltonian matrix and $J_{\mu\nu}=\sum_{\lambda,\sigma}P_{\lambda\sigma}\langle\mu\lambda|\nu\sigma\rangle$ is the Coulomb matrix, matching the standard Hartree–Fock operators. For hybrid density functionals, a fraction $a_{\text{x}}$ of Hartree–Fock exact exchange $K_{\mu\nu}=\sum_{\lambda,\sigma}P_{\lambda\sigma}\langle\mu\lambda|\sigma\nu\rangle$ is incorporated with spin-scaling factor $c_{\text{x}}$, while $V_{\mu\nu}^{\text{xc}}=\langle\chi_\mu|v_{\text{xc}}|\chi_\nu\rangle$ represents the matrix elements of the exchange-correlation potential.

The one-particle density matrix is constructed from the occupied orbital coefficients with orbital occupancy factor $f_{\text{occ}}$ as

$$
P_{\mu\nu}=f_{\text{occ}}\sum_i^{\text{occ}}C_{\mu i}C_{\nu i}^*
$$

where $f_{\text{occ}}=2$ and $c_{\text{x}}=1/2$ for Restricted Kohn–Sham (RKS) describing closed-shell systems with doubly occupied spatial orbitals, while $f_{\text{occ}}=1$ and $c_{\text{x}}=1$ for Generalized Kohn–Sham (GKS) where each state is an individually occupied two-component spinor. The total electronic energy of the molecular system is evaluated from the density matrix and potential contributions as

$$
E_{\text{elec}}=\sum_{\mu,\nu}P_{\nu\mu}H_{\mu\nu}+\frac{1}{2}\sum_{\mu,\nu}P_{\nu\mu}\left(J_{\mu\nu}-c_{\text{x}}a_{\text{x}}K_{\mu\nu}\right)+E_{\text{xc}}[\rho]
$$

and the total molecular ground-state energy is obtained by adding the classical nuclear repulsion energy $V_{\text{nuc}}=\sum_{A<B}\frac{Z_A Z_B}{|\mathbf{R}_A-\mathbf{R}_B|}$ as $E_{\text{tot}}=E_{\text{elec}}+V_{\text{nuc}}$.

---

## III. Exchange-Correlation Approximations and Numerical Quadrature

### 5. Jacob's Ladder of Density Functional Approximations

Because the exact mathematical form of $E_{\text{xc}}[\rho]$ is unknown, practical DFT relies on systematic approximations organized along Perdew's conceptual Jacob's ladder according to the local physical ingredients included across successive rungs.

The Local Density Approximation (LDA) forms the first rung, modeling the exchange-correlation energy density at each point from the analytical and Quantum Monte Carlo parameterization of the homogeneous electron gas as

$$
E_{\text{xc}}^{\text{LDA}}[\rho]=\int\rho(\mathbf{r})\epsilon_{\text{xc}}^{\text{unif}}(\rho(\mathbf{r}))d\mathbf{r}
$$

The Generalized Gradient Approximation (GGA) forms the second rung, incorporating the density gradient norm squared $\gamma=|\nabla\rho|^2$ through the dimensionless reduced density gradient $s=\frac{|\nabla\rho|}{2(3\pi^2)^{1/3}\rho^{4/3}}$ to describe the spatial inhomogeneity of molecular electron distributions, as exemplified by PBE and BLYP.

The meta-GGA level forms the third rung, adding the non-interacting positive kinetic energy density $\tau(\mathbf{r})=\frac{1}{2}\sum_i^{\text{occ}}|\nabla\phi_i(\mathbf{r})|^2$ or the density Laplacian $\nabla^2\rho$ to distinguish between single-bond, metallic, and non-covalent overlapping regions, as implemented in SCAN and TPSS.

Hybrid functionals occupy the fourth rung, mixing a fraction $a_{\text{x}}$ of non-local Hartree–Fock exact exchange into the functional as $E_{\text{xc}}^{\text{hybrid}}=E_{\text{xc}}^{\text{DFA}}+a_{\text{x}}\left(E_{\text{x}}^{\text{HF}}-E_{\text{x}}^{\text{DFA}}\right)$, which partially cancels the self-interaction error and substantially improves reaction barrier heights, band gaps, and thermochemical atomization energies, as in B3LYP and PBE0.

### 6. Molecular Numerical Grid Quadrature

Because exchange-correlation energy densities $\epsilon_{\text{xc}}(\rho,\gamma,\tau)$ are non-linear functions of the density and its gradients, spatial integrals involving $E_{\text{xc}}$ and $V_{\mu\nu}^{\text{xc}}$ cannot be evaluated analytically over Gaussian basis sets. Instead, the continuous integrals are evaluated using numerical quadrature over a three-dimensional molecular grid.

The molecular integration space is partitioned into atom-centered contributions using Becke fuzzy cell weighting functions $w_A(\mathbf{r})$, mapping the continuous exchange-correlation integral into a discrete quadrature sum over grid points as

$$
E_{\text{xc}}\approx\sum_g w_g\rho(\mathbf{r}_g)\epsilon_{\text{xc}}(\mathbf{r}_g)
$$

where $\mathbf{r}_g$ are the Cartesian coordinates of grid point $g$, and $w_g$ is the quadrature weight combining radial integration steps and angular Lebedev weights on atomic spheres. At each grid point $\mathbf{r}_g$, the density and its gradient are evaluated directly from the atomic basis functions and the density matrix as

$$
\rho(\mathbf{r}_g)=\sum_{\mu,\nu}P_{\mu\nu}\chi_\mu(\mathbf{r}_g)\chi_\nu(\mathbf{r}_g),\quad\nabla\rho(\mathbf{r}_g)=2\sum_{\mu,\nu}P_{\mu\nu}\nabla\chi_\mu(\mathbf{r}_g)\chi_\nu(\mathbf{r}_g)
$$

Applying the chain rule to the quadrature expression with respect to variations in the density matrix elements yields the exchange-correlation potential matrix elements as

$$
V_{\mu\nu}^{\text{xc}}\approx\sum_g w_g\left[\frac{\partial(\rho\epsilon_{\text{xc}})}{\partial\rho}\chi_\mu\chi_\nu+2\frac{\partial(\rho\epsilon_{\text{xc}})}{\partial\gamma}\nabla\rho\cdot\left(\nabla\chi_\mu\chi_\nu+\chi_\mu\nabla\chi_\nu\right)+\frac{1}{2}\frac{\partial(\rho\epsilon_{\text{xc}})}{\partial\tau}\nabla\chi_\mu\cdot\nabla\chi_\nu\right]
$$

which is assembled on the numerical grid and added to the core and Coulomb matrices during each self-consistent field iteration until orbital coefficients and electron density achieve numerical convergence.
