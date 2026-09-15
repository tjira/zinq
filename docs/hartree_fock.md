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

where the diagonal blocks $\mathbf{F}^{\alpha\alpha}$ and $\mathbf{F}^{\beta\beta}$ describe the energy and interactions that preserve the spin of the electrons, while the off-diagonal blocks $\mathbf{F}^{\alpha\beta}$ and $\mathbf{F}^{\beta\alpha}$ describe spin-mixing exchange interactions that couple alpha and beta channels. The matrices $\mathbf{C}^{\alpha}$ and $\mathbf{C}^{\beta}$ describe the spatial distribution of the alpha and beta spin components of the molecular orbitals, allowing the system to model non-collinear spin arrangements.

---

## II. Energy and SCF Optimization

### 3. Fock Matrix Construction

The elements of the Fock matrix in the basis functions are constructed as

$$
F_{\mu\nu}=H_{\mu\nu}^{\text{core}}+\sum_{\lambda,\sigma}P_{\lambda\sigma}\left(\langle\mu\lambda|\nu\sigma\rangle-c_{\text{x}}\langle\mu\lambda|\sigma\nu\rangle\right)
$$

where $H_{\mu\nu}^{\text{core}}$ is the core Hamiltonian matrix containing the kinetic energy of the electrons and their electrostatic attraction to the nuclei, $P_{\lambda\sigma}$ represents the elements of the density matrix, $\langle\mu\lambda|\nu\sigma\rangle$ represents the two-electron Coulomb integrals describing the classical electrostatic repulsion between electron clouds, and $\langle\mu\lambda|\sigma\nu\rangle$ represents the two-electron exchange integrals. The exchange term is a purely quantum mechanical effect arising from the antisymmetry of the wavefunction, scaled by the exchange factor $c_{\text{x}}$. The density matrix $\mathbf{P}$ is computed from the occupied molecular orbital coefficients as

$$
P_{\lambda\sigma}=f_{\text{occ}}\sum_i^{\text{occ}}C_{\lambda i}C_{\sigma i}
$$

accounting for the occupied orbitals with an orbital occupancy factor $f_{\text{occ}}$. In Restricted Hartree–Fock (RHF), the equations are solved in the spatial basis where each spatial orbital is doubly occupied ($f_{\text{occ}}=2$), and the exchange factor is $c_{\text{x}}=1/2$ because exchange only occurs between electrons of identical spin. In Generalized Hartree–Fock (GHF), the equations are formulated in the combined spin-orbital basis where each spin-orbital has single occupancy ($f_{\text{occ}}=1$), and the exchange factor is $c_{\text{x}}=1$ because spin is explicitly resolved within the basis.

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

where the factor of $1/2$ is necessary to prevent double-counting the electron-electron repulsions that are included in the Fock matrix. The total molecular energy is obtained by adding the classical electrostatic repulsion energy between the nuclei as

$$
E_{\text{tot}}=E_{\text{elec}}+V_{\text{nuc}}
$$

where the nuclear repulsion energy is evaluated from the atomic charges $Z_A$ and Cartesian nuclear positions $\mathbf{R}_A$ as

$$
V_{\text{nuc}}=\sum_{A<B}\frac{Z_AZ_B}{|\mathbf{R}_A-\mathbf{R}_B|}
$$

which completes the total ground-state energy evaluation.

---

## III. Nuclear Derivatives and Response

### 6. Analytical Nuclear Gradient

To find the forces acting on the atoms (which we need for moving atoms in molecular dynamics or finding stable geometries), we calculate the derivative of the Hartree–Fock energy with respect to the nuclear coordinates. Because the Hartree–Fock wavefunction is variationally optimized, its energy is stationary with respect to changes in the molecular orbital coefficients, meaning their derivatives do not appear in the gradient. The term containing the derivative of the overlap matrix $\mathbf{S}$ represents the Pulay force, which arises because the basis functions are centered on the atoms and move along with them as the nuclei move. The analytical nuclear gradient is given by

$$
\frac{dE_{\text{HF}}}{dx}=\frac{dV_{\text{nuc}}}{dx}+\sum_{\mu,\nu}P_{\mu\nu}\frac{dH_{\mu\nu}^{\text{core}}}{dx}-\sum_{\mu,\nu}W_{\mu\nu}\frac{dS_{\mu\nu}}{dx}+\frac{1}{2}\sum_{\mu,\nu,\lambda,\sigma}P_{\mu\nu}P_{\lambda\sigma}\left(\frac{d\langle\mu\lambda|\nu\sigma\rangle}{dx}-c_{\text{x}}\frac{d\langle\mu\lambda|\sigma\nu\rangle}{dx}\right)
$$

where $V_{\text{nuc}}$ is the nuclear repulsion energy, $c_{\text{x}}$ is the exchange scaling factor, and $\mathbf{W}$ is the energy-weighted density matrix defined from the molecular orbital energies $\epsilon_i$ and coefficients as

$$
W_{\mu\nu}=f_{\text{occ}}\sum_i^{\text{occ}}\epsilon_iC_{\mu i}C_{\nu i}
$$

accounting for the energy of the occupied orbitals with an orbital occupancy factor $f_{\text{occ}}$. In Restricted Hartree–Fock (RHF), the calculation is performed in the spatial atomic orbital basis where each spatial orbital is doubly occupied ($f_{\text{occ}}=2$), and the exchange factor is $c_{\text{x}}=1/2$ because exchange only occurs between electrons of identical spin. In Generalized Hartree–Fock (GHF), the calculation is formulated in a combined spin-orbital basis of twice the spatial dimension where each spin-orbital has single occupancy ($f_{\text{occ}}=1$), and the exchange factor is $c_{\text{x}}=1$ because the spin integration is carried out directly over the spin-orbitals.

### 7. Coupled-Perturbed Hartree–Fock Equations

When a molecular geometry is perturbed, the molecular orbitals rotate to preserve the self-consistent field condition $\mathbf{F}\mathbf{C}=\mathbf{S}\mathbf{C}\mathbf{E}$. The Coupled-Perturbed Hartree–Fock (CPHF) equations solve directly for the orbital response matrix $\mathbf{U}^x$, which describes the first-order transformation of the molecular orbital coefficients as

$$
\frac{dC_{\mu p}}{dx}=\sum_qC_{\mu q}U_{qp}^x
$$

where the occupied–virtual blocks $U_{ai}^x$ represent the physical relaxation and polarization of the electronic wavefunction, while the remaining blocks preserve orbital orthonormality through the constraint

$$
U_{pq}^x+U_{qp}^x+S_{pq}^{x,\text{MO}}=0
$$

where $S_{pq}^{x,\text{MO}}$ denotes the overlap derivative in the molecular orbital basis. Because rotating the orbitals induces a first-order perturbed density matrix $\mathbf{P}^x$ that modifies the Fock operator through electron repulsion, the response matrix elements $U_{ai}^x$ are coupled and solved iteratively as

$$
(\epsilon_a-\epsilon_i)U_{ai}^x=-\left(F_{ai}^{x,\text{MO}}+V_{ai}^{x,\text{MO}}(\mathbf{P}^x)-S_{ai}^{x,\text{MO}}\epsilon_i\right)
$$

where $F_{ai}^{x,\text{MO}}$ is the MO-transformed skeleton Fock derivative, and $V_{ai}^{x,\text{MO}}(\mathbf{P}^x)$ is the response potential generated by the first-order perturbed density matrix

$$
P_{\mu\nu}^x=f_{\text{occ}}\sum_j^{\text{occ}}\left(\frac{dC_{\mu j}}{dx}C_{\nu j}+C_{\mu j}\frac{dC_{\nu j}}{dx}\right)
$$

which is contracted with the two-electron integrals using an exchange scaling factor $c_{\text{x}}$ as

$$
V_{\lambda\sigma}^x=\sum_{\mu,\nu}P_{\mu\nu}^x\left(\langle\mu\lambda|\nu\sigma\rangle-c_{\text{x}}\langle\mu\nu|\lambda\sigma\rangle\right)
$$

before transformation into the molecular orbital basis. This formulation is universal for both spin variants, where Restricted Hartree–Fock uses $f_{\text{occ}}=2$ and $c_{\text{x}}=1/2$ in the spatial basis, whereas Generalized Hartree–Fock uses $f_{\text{occ}}=1$ and $c_{\text{x}}=1$ in the spin-orbital basis. Once the linear equations are solved for the response matrix $\mathbf{U}^x$, the molecular orbital coefficient derivatives are obtained by back-transformation as $\frac{d\mathbf{C}}{dx}=\mathbf{C}\mathbf{U}^x$, and the orbital energy derivatives are evaluated from the diagonal elements of the perturbed Fock matrix as

$$
\frac{d\epsilon_p}{dx}=F_{pp}^{x,\text{MO}}+V_{pp}^{x,\text{MO}}(\mathbf{P}^x)-S_{pp}^{x,\text{MO}}\epsilon_p
$$

which describes the first-order shift in each orbital energy eigenvalue. While first derivatives of the variational Hartree–Fock energy do not require orbital response calculations due to the Hellmann–Feynman theorem, the solved response matrix $\mathbf{U}^x$ and orbital derivatives are essential for calculating second derivatives of the Hartree–Fock energy (nuclear Hessians and vibrational frequencies) as well as analytical nuclear gradients for non-variational correlated methods such as Møller–Plesset perturbation theory.
