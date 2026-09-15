# Density Functional Theory

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of Density Functional Theory (DFT). If you want to understand how the electronic structure of molecules can be determined from the three-dimensional electron density rather than a many-electron wavefunction, this guide is written step-by-step for you.

---

## I. Kohn–Sham Formulation

In quantum mechanics, the exact wavefunction for $N$ electrons depends simultaneously on $3N$ spatial coordinates, making exact solutions computationally intractable for larger systems. Density Functional Theory resolves this challenge through the Hohenberg–Kohn theorems, which prove that the ground-state properties and total energy of an electronic system are uniquely determined by the three-dimensional electronic density $\rho(\mathbf{r})$.

### 1. Kohn–Sham Equations

The Kohn–Sham formulation maps the interacting electron system onto an auxiliary system of non-interacting electrons experiencing an effective local potential, yielding the Roothaan–Kohn–Sham matrix equation

$$
\mathbf{F}\mathbf{C}=\mathbf{S}\mathbf{C}\mathbf{E}
$$

where $\mathbf{F}$ is the Kohn–Sham Fock matrix, $\mathbf{C}$ is the molecular orbital coefficient matrix, $\mathbf{S}$ is the overlap matrix of the atomic orbital basis functions, and $\mathbf{E}$ is the diagonal matrix containing the orbital energy eigenvalues. In the atomic orbital basis, the Kohn–Sham matrix elements are constructed as

$$
F_{\mu\nu}=H_{\mu\nu}^{\text{core}}+\sum_{\lambda,\sigma}P_{\lambda\sigma}\left(\langle\mu\lambda|\nu\sigma\rangle-c_{\text{x}}a_{\text{x}}\langle\mu\lambda|\sigma\nu\rangle\right)+V_{\mu\nu}^{\text{xc}}
$$

where $H_{\mu\nu}^{\text{core}}$ is the one-electron core Hamiltonian matrix describing kinetic energy and nuclear attraction, the two-electron integral sum accounts for classical Coulomb repulsion and an exact exchange contribution scaled by the hybrid mixing fraction $a_{\text{x}}$ and spin-variant exchange factor $c_{\text{x}}$, and $V_{\mu\nu}^{\text{xc}}$ is the exchange-correlation potential matrix. The density matrix $\mathbf{P}$ is formed from the occupied molecular orbital coefficients, while the exchange-correlation potential encapsulates all non-classical electron-electron interactions.

### 2. Restricted and Generalized Variants

Depending on how electron spin is treated, Kohn–Sham density functional theory is formulated in either a restricted or generalized framework. The density matrix is generally constructed from the occupied molecular orbital coefficients with an orbital occupancy factor $f_{\text{occ}}$ as

$$
P_{\mu\nu}=f_{\text{occ}}\sum_i^{\text{occ}} C_{\mu i}C_{\nu i}^*
$$

where the occupancy factor reflects the nature of the basis. In Restricted Kohn–Sham (RKS), electrons of opposite spins are constrained to share identical spatial molecular orbitals, yielding a spin-unpolarized density $\rho_\alpha(\mathbf{r})=\rho_\beta(\mathbf{r})=\frac{1}{2}\rho(\mathbf{r})$. The density matrix $\mathbf{P}$ is defined in the $N_{\text{ao}}$-dimensional spatial basis where each spatial orbital holds two electrons ($f_{\text{occ}}=2$) with exchange factor $c_{\text{x}}=1/2$. In Generalized Kohn–Sham (GKS), the molecular orbitals are expanded into a $2N_{\text{ao}}$-dimensional spin-orbital basis where each state is a general two-component spinor with single occupancy ($f_{\text{occ}}=1$) and exchange factor $c_{\text{x}}=1$. The generalized density matrix is partitioned into spin blocks as

$$
\mathbf{P}=\begin{pmatrix}\mathbf{P}^{\alpha\alpha} & \mathbf{P}^{\alpha\beta} \\ \mathbf{P}^{\beta\alpha} & \mathbf{P}^{\beta\beta}\end{pmatrix}
$$

where the diagonal blocks $\mathbf{P}^{\alpha\alpha}$ and $\mathbf{P}^{\beta\beta}$ determine the individual spin densities $\rho_\alpha(\mathbf{r})$ and $\rho_\beta(\mathbf{r})$ on the numerical grid, while the off-diagonal blocks $\mathbf{P}^{\alpha\beta}$ and $\mathbf{P}^{\beta\alpha}$ describe non-collinear spin alignments that couple through the exact exchange operator in hybrid functionals. The resulting exchange-correlation potential matrix $\mathbf{V}^{\text{xc}}$ is block-diagonal in the spin channels with diagonal blocks $\mathbf{V}^{\text{xc},\alpha}$ and $\mathbf{V}^{\text{xc},\beta}$ evaluated from the spin-polarized functional derivatives.

### 3. Exchange-Correlation and Total Energy

The exchange-correlation energy accounts for the quantum mechanical exchange and dynamical correlation effects, evaluated as the spatial integral of the exchange-correlation energy density $\epsilon_{\text{xc}}$ multiplied by the electronic density as

$$
E_{\text{xc}}=\int\rho(\mathbf{r})\epsilon_{\text{xc}}(\mathbf{r})d\mathbf{r}
$$

where the functional form of $\epsilon_{\text{xc}}$ depends on the chosen density functional approximation. The total electronic energy is then calculated from the density matrix, core Hamiltonian, Kohn–Sham matrix, and exchange-correlation potential as

$$
E_{\text{elec}}=\frac{1}{2}\sum_{\mu,\nu}P_{\mu\nu}\left(H_{\mu\nu}^{\text{core}}+F_{\mu\nu}-V_{\mu\nu}^{\text{xc}}\right)+E_{\text{xc}}
$$

and the total energy of the molecular system is obtained by adding the classical nuclear repulsion energy $V_{\text{nuc}}$ as

$$
E_{\text{tot}}=E_{\text{elec}}+V_{\text{nuc}}
$$

where $V_{\text{nuc}}$ is evaluated from the nuclear charges and Cartesian nuclear coordinates.

---

## II. Numerical Grid Integration

### 4. Molecular Integration Grid

Because the exchange-correlation energy density $\epsilon_{\text{xc}}$ is a highly non-linear function of the density and its spatial derivatives, the spatial integrals cannot be evaluated analytically. Instead, the integration is performed numerically over a three-dimensional grid of points distributed across the molecule. The molecular volume is partitioned into atomic regions using Becke fuzzy cell weighting, where a multicenter grid is formed for each atom by combining a radial grid with an angular Lebedev quadrature sphere. The exchange-correlation energy is evaluated as a weighted quadrature sum over all grid points as

$$
E_{\text{xc}}\approx\sum_gw_g\rho(\mathbf{r}_g)\epsilon_{\text{xc}}(\mathbf{r}_g)
$$

where $\mathbf{r}_g$ represents the spatial coordinates of grid point $g$, and $w_g$ is the quadrature weight incorporating the radial, angular, and multicenter Becke partitioning factors.

### 5. Density Evaluation and Functional Families

At each quadrature grid point, the electronic density $\rho(\mathbf{r})$ is evaluated from the atomic orbital basis functions $\chi_\mu(\mathbf{r})$ and density matrix $\mathbf{P}$ as

$$
\rho(\mathbf{r})=\sum_{\mu,\nu}P_{\mu\nu}\chi_\mu(\mathbf{r})\chi_\nu(\mathbf{r})
$$

along with the density gradient $\nabla\rho(\mathbf{r})$ and the non-interacting kinetic energy density $\tau(\mathbf{r})$ evaluated as

$$
\nabla\rho(\mathbf{r})=2\sum_{\mu,\nu}P_{\mu\nu}\nabla\chi_\mu(\mathbf{r})\chi_\nu(\mathbf{r})
$$

and

$$
\tau(\mathbf{r})=\frac{1}{2}\sum_{\mu,\nu}P_{\mu\nu}\nabla\chi_\mu(\mathbf{r})\cdot\nabla\chi_\nu(\mathbf{r})
$$

which classify functionals into standard approximations: the Local Density Approximation (LDA) depends solely on the local density $\rho$; the Generalized Gradient Approximation (GGA) incorporates the gradient norm squared $\gamma=|\nabla\rho|^2$; meta-GGA functionals additionally include the kinetic energy density $\tau$ or density Laplacian $\nabla^2\rho$; and hybrid functionals mix a fraction $a_{\text{x}}$ of non-local Hartree–Fock exact exchange into the Kohn–Sham operator. In spin-polarized generalized calculations, these expressions are evaluated separately for the alpha and beta spin channels using the diagonal blocks $\mathbf{P}^{\alpha\alpha}$ and $\mathbf{P}^{\beta\beta}$.

### 6. Exchange-Correlation Potential Construction

The matrix elements of the exchange-correlation potential $V_{\mu\nu}^{\text{xc}}$ are defined as the functional derivative of the exchange-correlation energy with respect to the density matrix elements. Applying the chain rule to the numerical quadrature expansion over grid points yields the matrix representation as

$$
V_{\mu\nu}^{\text{xc}}\approx\sum_gw_g\left(\frac{\partial(\rho\epsilon_{\text{xc}})}{\partial\rho}\chi_\mu(\mathbf{r}_g)\chi_\nu(\mathbf{r}_g)+2\frac{\partial(\rho\epsilon_{\text{xc}})}{\partial\gamma}\nabla\rho(\mathbf{r}_g)\cdot\left(\nabla\chi_\mu(\mathbf{r}_g)\chi_\nu(\mathbf{r}_g)+\chi_\mu(\mathbf{r}_g)\nabla\chi_\nu(\mathbf{r}_g)\right)+\frac{1}{2}\frac{\partial(\rho\epsilon_{\text{xc}})}{\partial\tau}\nabla\chi_\mu(\mathbf{r}_g)\cdot\nabla\chi_\nu(\mathbf{r}_g)\right)
$$

where the partial derivatives of the exchange-correlation energy density are contracted with the basis function values and their spatial gradients at each grid point to assemble the potential matrix for self-consistent field iterations.
