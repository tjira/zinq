# Density Functional Theory

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of Density Functional Theory (DFT) and how it is implemented in our framework. If you have ever wondered how scientists can simulate the quantum mechanical behavior of hundreds of electrons in a molecule without running out of computer memory, this guide is written step-by-step for you.

---

## I. Kohn–Sham Formulation

To understand Density Functional Theory, we must first understand why simulating many electrons is hard. In quantum mechanics, if we have $N$ electrons, the wavefunction that describes them depends on all their coordinates simultaneously. This means that as we add more electrons, the mathematical complexity of the problem grows exponentially, making it impossible to solve for larger molecules. Density Functional Theory solves this problem by using a simple but powerful idea: instead of trying to find a high-dimensional wavefunction, we can find the 3D electronic density, which represents the probability of finding any electron at a given point in space. The Hohenberg–Kohn theorems prove that the ground-state properties of a system are completely determined by this 3D electronic density, which dramatically simplifies the calculation.

### 1. Kohn–Sham Equations

The Kohn–Sham formulation of Density Functional Theory maps the complicated system of interacting electrons onto an auxiliary, simplified system of non-interacting electrons that move through an effective local potential. This formulation leads to a set of matrix equations called the Kohn–Sham equations

$$
\mathbf{F}\mathbf{C}=\mathbf{S}\mathbf{C}\mathbf{E}
$$

where $\mathbf{F}$ is the Kohn–Sham matrix (which acts like the Fock matrix in Hartree–Fock theory), $\mathbf{C}$ is the molecular orbital coefficient matrix, $\mathbf{S}$ is the overlap matrix of the atomic orbital basis functions, and $\mathbf{E}$ is the diagonal matrix containing the molecular orbital energy eigenvalues. In the atomic orbital basis, the elements of the Kohn–Sham matrix are constructed as

$$
F_{\mu\nu}=H_{\mu\nu}^{\text{core}}+J_{\mu\nu}+V_{\mu\nu}^{\text{xc}}
$$

where $H_{\mu\nu}^{\text{core}}$ is the core Hamiltonian matrix containing the kinetic energy of the electrons and their electrostatic attraction to the atomic nuclei, $J_{\mu\nu}$ is the classical Coulomb repulsion matrix describing how the electronic density repels itself, and $V_{\mu\nu}^{\text{xc}}$ is the exchange-correlation potential matrix. The exchange-correlation potential is the most important term in DFT because it contains all the quantum mechanical effects, such as exchange (the Pauli exclusion principle that keeps electrons of the same spin apart) and correlation (the detailed movements of electrons to avoid each other), that are not captured by the classical terms.

### 2. Exchange-Correlation Potential and Energy

The total exchange-correlation energy is calculated by integrating the exchange-correlation energy density per unit volume $\epsilon_{\text{xc}}$ multiplied by the electronic density $\rho$ over all spatial coordinates as

$$
E_{\text{xc}}=\int\rho(\mathbf{r})\epsilon_{\text{xc}}(\rho(\mathbf{r}))d\mathbf{r}
$$

where $\mathbf{r}$ is the 3D position vector and $\rho(\mathbf{r})$ is the electronic density at that position. The exchange-correlation potential matrix elements $V_{\mu\nu}^{\text{xc}}$ are then obtained by taking the functional derivative of this exchange-correlation energy with respect to the electronic density.

---

## II. Numerical Grid Integration

### 3. Molecular Grid Construction

Because the exchange-correlation energy density $\epsilon_{\text{xc}}$ is a highly non-linear and mathematically complex function of the density, it is impossible to solve the integration analytically with pen and paper. Instead, we must perform the integration numerically over a three-dimensional grid of points surrounding the molecule. To construct this grid, the program partitions the space around the molecule into individual atomic regions using Becke's fuzzy cell partitioning scheme. For each atom, a spherical grid is constructed by combining a radial grid (which extends outwards from the nucleus using Euler–Maclaurin or Treutler quadrature schemes) and an angular grid (which distributes points on a sphere around the nucleus using Lebedev quadrature). The total integration is then evaluated as a weighted sum over all grid points as

$$
E_{\text{xc}}\approx\sum_gw_g\rho(\mathbf{r}_g)\epsilon_{\text{xc}}(\rho(\mathbf{r}_g))
$$

where $\mathbf{r}_g$ represents the spatial coordinates of grid point $g$ and $w_g$ is the corresponding Becke-partitioned integration weight that represents the volume element associated with that grid point.

### 4. Density Evaluation and Functional Families

At each grid point, the electronic density and its spatial derivatives are evaluated by contracting the density matrix $\mathbf{P}$ with the atomic orbital basis functions $\chi_\mu$ and their derivatives. The exchange-correlation energy and potential densities are then evaluated using the libxc library. Functionals in DFT are grouped into different families based on the variables they depend on: the Local Density Approximation (LDA) family depends only on the local electronic density $\rho$; the Generalized Gradient Approximation (GGA) family also depends on the local gradient of the density, represented by the gradient norm squared $\gamma=|\nabla\rho|^2$, which accounts for how fast the density changes in space; and the meta-GGA family further incorporates the local kinetic energy density $\tau$, which represents the kinetic energy of the non-interacting reference system and allows for highly accurate calculations.

### 5. Implementation Tricks and Memory Structures

Evaluating the basis functions and density at millions of grid points is computationally expensive. To optimize this process, our codebase evaluates the basis function values and their spatial derivatives at all grid points once at the start of the calculation and stores them in a memory structure named `BasisGrid(T)`. During each self-consistent field iteration, the electronic density and its gradients are computed at each grid point by contracting the density matrix $\mathbf{P}$ with the pre-evaluated basis functions, storing the results in a structure named `DensityGrid(T)`. The resulting energy and potential values calculated by libxc are collected in a structure named `PotentialGrid(T)` and contracted back to the atomic orbital basis to form the exchange-correlation potential matrix $\mathbf{V}_{\text{xc}}$. This structured approach avoids re-evaluating basis functions on the grid points during each iteration, significantly accelerating the simulation.
