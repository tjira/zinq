# Hartree–Fock Method

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of the Hartree–Fock method and how it is implemented in our framework. If you want to understand the fundamental physics and algorithms that form the starting point of modern quantum chemistry calculations, this guide is written step-by-step for you.

---

## I. Core Wavefunction Theory

To understand the Hartree–Fock method, we must first understand how it simplifies the description of electrons in a molecule. In quantum mechanics, electrons are identical particles that carry spin and repel each other. The behavior of these electrons is described by a many-electron wavefunction. Because electrons are fermions, their wavefunction must be antisymmetric, which means that if we swap the positions of any two electrons, the wavefunction must change its sign. The Hartree–Fock method models this many-electron wavefunction as a single mathematical construct called a Slater determinant. A Slater determinant is built from individual one-electron functions called molecular orbitals. By using a determinant, the wavefunction automatically satisfies both the Pauli exclusion principle (no two electrons can occupy the exact same state) and the antisymmetry requirement. In this model, we make a mean-field approximation: we assume that each electron does not interact with the other electrons individually, but instead moves through an average electric field created by all the other electrons combined.

### 1. Roothaan–Hall Equations

To solve the Hartree–Fock equations on a computer, we express the molecular orbitals as a linear combination of simpler, known functions centered on the atoms, which we call atomic orbital basis functions. When we apply the variational principle to find the set of molecular orbitals that minimizes the total electronic energy, we obtain a set of matrix equations called the Roothaan–Hall equations

$$
\mathbf{F}\mathbf{C}=\mathbf{S}\mathbf{C}\mathbf{E}
$$

where $\mathbf{F}$ is the Fock matrix, which represents the effective energy operator for a single electron moving in the average field of all others, $\mathbf{C}$ is the molecular orbital coefficient matrix that tells us how to combine the atomic orbitals to form the molecular orbitals, $\mathbf{S}$ is the overlap matrix that measures how much the non-orthogonal atomic orbital basis functions overlap with each other in space, and $\mathbf{E}$ is a diagonal matrix containing the molecular orbital energies. Because the atomic basis functions are not orthogonal to each other, the overlap matrix $\mathbf{S}$ is not the identity matrix. This turns the problem into a generalized eigenvalue problem rather than a standard one. To solve these equations, we must first transform the matrices into an orthogonal basis where the overlap matrix becomes the identity matrix, solve the standard eigenvalue problem, and then transform the resulting coefficients back. In the generalized variant of Hartree–Fock, we use spin-orbitals that explicitly include the spin state (alpha or beta) of the electron, which doubles the size of our matrices.

### 2. Generalized Hartree–Fock and Spin Variants

Depending on how we handle the spin of the electrons, there are different versions of the Hartree–Fock method: Restricted Hartree–Fock (RHF) forces electrons of opposite spins (alpha and beta) to share the exact same spatial molecular orbitals, which works well for stable molecules with paired electrons but fails when we try to break chemical bonds; Unrestricted Hartree–Fock (UHF) allows alpha and beta electrons to have different spatial orbitals; and Generalized Hartree–Fock (GHF) goes even further by allowing each molecular orbital to be a general mixture of alpha and beta spin components, meaning the spin direction can vary continuously in space. The GHF equations are formulated in a combined spin-orbital basis of twice the spatial dimension as

$$
\begin{pmatrix}\mathbf{F}^{\alpha\alpha}&\mathbf{F}^{\alpha\beta}\\\mathbf{F}^{\beta\alpha}&\mathbf{F}^{\beta\beta}\end{pmatrix}\begin{pmatrix}\mathbf{C}^{\alpha}\\\mathbf{C}^{\beta}\end{pmatrix}=\begin{pmatrix}\mathbf{S}&\mathbf{0}\\\mathbf{0}&\mathbf{S}\end{pmatrix}\begin{pmatrix}\mathbf{C}^{\alpha}\\\mathbf{C}^{\beta}\end{pmatrix}\mathbf{E}
$$

where the diagonal blocks $\mathbf{F}^{\alpha\alpha}$ and $\mathbf{F}^{\beta\beta}$ describe the energy and interactions that preserve the spin of the electrons, while the off-diagonal blocks $\mathbf{F}^{\alpha\beta}$ and $\mathbf{F}^{\beta\alpha}$ describe spin-mixing interactions (such as spin-orbit coupling). The matrices $\mathbf{C}^{\alpha}$ and $\mathbf{C}^{\beta}$ describe the spatial distribution of the alpha and beta spin components of the molecular orbitals, allowing the system to model complex magnetic arrangements.

---

## II. Energy and SCF Optimization

### 3. Fock Matrix Construction

For a closed-shell system where all electrons are paired in spatial orbitals (Restricted Hartree–Fock), the elements of the Fock matrix in the atomic orbital basis are constructed as

$$
F_{\mu\nu}=H_{\mu\nu}^{\text{core}}+\sum_{\lambda,\sigma}P_{\lambda\sigma}\left(\langle\mu\lambda|\nu\sigma\rangle-\frac{1}{2}\langle\mu\lambda|\sigma\nu\rangle\right)
$$

where $H_{\mu\nu}^{\text{core}}$ is the core Hamiltonian matrix containing the kinetic energy of the electrons and their electrostatic attraction to the nuclei, $P_{\lambda\sigma}$ represents the elements of the density matrix, $\langle\mu\lambda|\nu\sigma\rangle$ represents the two-electron Coulomb integrals describing the classical electrostatic repulsion between electron clouds, and $\langle\mu\lambda|\sigma\nu\rangle$ represents the two-electron exchange integrals. The exchange term is a purely quantum mechanical effect arising from the antisymmetry of the wavefunction, which acts to keep electrons of the same spin apart. The factor of $1/2$ on the exchange term arises because the summation runs over spatial orbitals, each of which can hold two electrons of opposite spins. The density matrix $\mathbf{P}$ is computed from the occupied molecular orbital coefficients as

$$
P_{\lambda\sigma}=2\sum_i^{\text{occ}}C_{\lambda i}C_{\sigma i}
$$

where the factor of two accounts for the double occupancy of each spatial orbital, and the sum runs over all occupied molecular orbitals. In Generalized Hartree–Fock, the equations are written directly in terms of spin-orbitals, which removes the factors of two and the $1/2$ scale factor.

### 4. Self-Consistent Field Iteration and DIIS

Because the Fock matrix $\mathbf{F}$ depends on the density matrix $\mathbf{P}$ (which is built from the molecular orbitals, which are the eigenvectors of $\mathbf{F}$), the Roothaan–Hall equations are non-linear and cannot be solved in a single step. Instead, we must solve them iteratively using the Self-Consistent Field (SCF) method. We start with a guess for the density matrix, use it to build a Fock matrix, diagonalize the Fock matrix to find new molecular orbitals, use these new orbitals to build a new density matrix, and repeat this cycle until the density matrix stops changing. To accelerate this process and prevent the calculation from oscillating or failing to converge, our codebase uses the Direct Inversion in the Iterative Subspace (DIIS) method. DIIS monitors the convergence by calculating an error matrix defined as

$$
\mathbf{e}=\mathbf{F}\mathbf{P}\mathbf{S}-\mathbf{S}\mathbf{P}\mathbf{F}
$$

where the error matrix $\mathbf{e}$ must become exactly zero at convergence, representing that the Fock matrix and density matrix commute in the orthogonalized basis. DIIS collects a history of Fock and error matrices from previous steps and solves a small system of linear equations with Lagrange multipliers to find the optimal weights to combine them into an accelerated Fock matrix for the next step.

### 5. Total Energy

Once the SCF cycle has converged, the total electronic Hartree–Fock energy of the system is calculated as

$$
E_{\text{elec}}=\frac{1}{2}\sum_{\mu,\nu}P_{\mu\nu}\left(H_{\mu\nu}^{\text{core}}+F_{\mu\nu}\right)
$$

where the factor of $1/2$ is necessary to prevent double-counting the electron-electron repulsions that are included in the Fock matrix. The total energy of the molecule is then the sum of this electronic energy and the classical electrostatic repulsion energy between the nuclei.

---

## III. Nuclear Derivatives and Response

### 6. Analytical Nuclear Gradient

To find the forces acting on the atoms (which we need for moving atoms in molecular dynamics or finding stable geometries), we calculate the derivative of the Hartree–Fock energy with respect to the nuclear coordinates. The analytical nuclear gradient is given by

$$
\frac{dE_{\text{HF}}}{dx}=\frac{dV_{\text{nuc}}}{dx}+\sum_{\mu,\nu}P_{\mu\nu}\frac{dH_{\mu\nu}^{\text{core}}}{dx}-\sum_{\mu,\nu}W_{\mu\nu}\frac{dS_{\mu\nu}}{dx}+\frac{1}{2}\sum_{\mu,\nu,\lambda,\sigma}P_{\mu\nu}P_{\lambda\sigma}\left(\frac{d\langle\mu\lambda|\nu\sigma\rangle}{dx}-\frac{1}{2}\frac{d\langle\mu\lambda|\sigma\nu\rangle}{dx}\right)
$$

where $V_{\text{nuc}}$ is the nuclear repulsion energy, and $\mathbf{W}$ is the energy-weighted density matrix defined from the molecular orbital energies $\epsilon_i$ and coefficients as

$$
W_{\mu\nu}=2\sum_i^{\text{occ}}\epsilon_iC_{\mu i}C_{\nu i}
$$

which accounts for the energy of the occupied orbitals. Because the Hartree–Fock wavefunction is variationally optimized, its energy is stationary with respect to changes in the molecular orbital coefficients, meaning their derivatives do not appear in the gradient. The term containing the derivative of the overlap matrix $\mathbf{S}$ represents the Pulay force, which arises because the basis functions are centered on the atoms and move along with them as the nuclei move.

### 7. Coupled-Perturbed Hartree–Fock Equations

When a molecule is perturbed, such as when an atom moves or when we apply an external electric field, the molecular orbitals change. Because the Fock matrix depends on the density matrix, any change in the orbitals affects the Fock matrix, which in turn affects the orbitals. The Coupled-Perturbed Hartree–Fock (CPHF) equations describe this coupled response self-consistently. We express the derivative of the molecular orbital coefficients in terms of the unperturbed coefficients using an orbital response matrix $\mathbf{U}^x$ as

$$
\frac{dC_{\mu i}}{dx}=\sum_pC_{\mu p}U_{pi}^x
$$

where the sum runs over all occupied and virtual molecular orbitals. The requirement that the molecular orbitals remain orthonormal as they change constrains the symmetric part of the response matrix to satisfy

$$
U_{pq}^x+U_{qp}^x+\frac{dS_{pq}}{dx}=0
$$

where $S_{pq}$ is the derivative of the overlap matrix in the molecular orbital basis. The remaining occupied-virtual blocks of the response matrix are found by solving the CPHF equations

$$
(\epsilon_a-\epsilon_i)U_{ai}^x-\sum_{j}^{\text{occ}}\sum_{b}^{\text{vir}}A_{ai,bj}U_{bj}^x=B_{ai}^x
$$

where $i, j$ denote occupied orbitals, $a, b$ denote virtual (unoccupied) orbitals, $A_{ai,bj}$ are the coupling matrix elements that describe how the Hartree–Fock potential changes when the density matrix is modified, and $B_{ai}^x$ is the direct perturbation vector containing the derivatives of the Fock and overlap matrices. Because building and inverting the coupling matrix $\mathbf{A}$ directly would require huge amounts of memory and time, our codebase solves these equations iteratively using DIIS acceleration.
