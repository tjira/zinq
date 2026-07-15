# Density Functional Theory

Density Functional Theory (DFT) is a quantum mechanical modeling method used to investigate the electronic structure of many-body systems. Instead of solving the many-electron Schrödinger equation for the wavefunction, DFT determines the properties of a system using the electronic density, which simplifies the computational complexity.

---

## I. Kohn–Sham Formulation

### 1. Kohn–Sham Equations

The Kohn–Sham formulation of density functional theory maps the interacting many-electron system onto an auxiliary system of non-interacting electrons moving in an effective local potential. This leads to the Kohn–Sham equations

$$
\mathbf{F}\mathbf{C}=\mathbf{S}\mathbf{C}\mathbf{E}
$$

where the Fock matrix $\mathbf{F}$ is replaced by the Kohn–Sham matrix. In the atomic orbital basis, the Kohn–Sham matrix is constructed as

$$
F_{\mu\nu}=H_{\mu\nu}^{\text{core}}+J_{\mu\nu}+V_{\mu\nu}^{\text{xc}}
$$

where $H_{\mu\nu}^{\text{core}}$ is the core Hamiltonian, $J_{\mu\nu}$ is the classical Coulomb potential of the electron density, and $V_{\mu\nu}^{\text{xc}}$ is the exchange-correlation potential matrix representing quantum mechanical exchange and correlation effects.

### 2. Exchange-Correlation Potential and Energy

The exchange-correlation energy is calculated by integrating the exchange-correlation energy density $\epsilon_{\text{xc}}$ over the spatial electronic density $\rho$ as

$$
E_{\text{xc}}=\int\rho(\mathbf{r})\epsilon_{\text{xc}}(\rho(\mathbf{r}))d\mathbf{r}
$$

where the exchange-correlation potential $V_{\mu\nu}^{\text{xc}}$ is obtained by taking the functional derivative of the exchange-correlation energy with respect to the density.

---

## II. Numerical Grid Integration

### 3. Molecular Grid Construction

Because the exchange-correlation energy density is a complex functional of the density, the integration is performed numerically over a molecular grid. The molecular grid is constructed by partitioning the molecular space into atomic regions using Becke's fuzzy cell partitioning scheme. For each atom, the grid is built by combining a radial grid (using Euler–Maclaurin or Treutler quadrature) and an angular grid (using Lebedev quadrature). The total integration is evaluated as a weighted sum over all grid points as

$$
E_{\text{xc}}\approx\sum_gw_g\rho(\mathbf{r}_g)\epsilon_{\text{xc}}(\rho(\mathbf{r}_g))
$$

where $\mathbf{r}_g$ represents the spatial coordinates of grid point $g$ and $w_g$ is the corresponding Becke-partitioned integration weight.

### 4. Density Evaluation and Functional Families

At each grid point, the electronic density and its derivatives are evaluated using the density matrix $\mathbf{P}$ and the atomic orbital basis functions $\chi_\mu$. The exchange-correlation potential and energy densities are evaluated using the libxc library, which classifies functionals into different families based on their variables. The Local Density Approximation (LDA) family depends only on the local electronic density, whereas the Generalized Gradient Approximation (GGA) family also incorporates the density gradient norm squared $\gamma=|\nabla\rho|^2$. The meta-GGA family further includes the kinetic energy density $\tau$ representing the local kinetic energy contribution, allowing highly accurate semi-local functional approximations to be constructed.

### 5. Implementation Tricks and Memory Structures

To optimize grid evaluations during the self-consistent field iterations, the codebase evaluates basis function values and their spatial derivatives on the integration grid points once and stores them in a `BasisGrid(T)` structure. During each SCF cycle, the values of $\rho$, $\gamma$, and $\tau$ are computed at each grid point by contracting the density matrix $\mathbf{P}$ with the pre-evaluated basis functions in a `DensityGrid(T)` structure. The resulting functional outputs from libxc are collected in a `PotentialGrid(T)` structure and mapped back to the atomic orbital basis via a fast matrix contraction to form the exchange-correlation potential matrix $\mathbf{V}_{\text{xc}}$, which avoids recomputing basis function values at grid points in each iteration.
