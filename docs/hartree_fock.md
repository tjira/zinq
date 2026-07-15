# Hartree–Fock Method

The Hartree–Fock method is a fundamental variational approach in quantum chemistry and electronic structure theory to approximate the ground-state electronic wavefunction of a multi-electron system. It models the wavefunction as a single Slater determinant of molecular orbitals, effectively introducing a mean-field approximation where each electron moves in the average electrostatic field of all other electrons.

---

## I. Core Wavefunction Theory

### 1. Roothaan–Hall Equations

By expanding the molecular orbitals as a linear combination of atomic orbitals, the variational optimization of the energy leads to the Roothaan–Hall equations

$$
\mathbf{F}\mathbf{C}=\mathbf{S}\mathbf{C}\mathbf{E}
$$

where $\mathbf{F}$ is the Fock matrix representing the effective one-electron Hamiltonian, $\mathbf{C}$ is the molecular orbital coefficient matrix relating the molecular orbitals to the atomic orbital basis functions, $\mathbf{S}$ is the atomic orbital overlap matrix representing the non-orthogonality of the basis functions, and $\mathbf{E}$ is the diagonal matrix of molecular orbital energies. The presence of the overlap matrix $\mathbf{S}$ defines this relation as a generalized eigenvalue problem rather than a standard one, arising directly from the use of a non-orthogonal atomic orbital basis set. For the generalized variant, these equations are solved using spin-orbitals instead of spatial orbitals, which doubles the dimension of the matrices to explicitly accommodate both alpha and beta spin components.

### 2. Generalized Hartree–Fock and Spin Variants

Generalized Hartree–Fock (GHF) generalizes the spin-orbital framework by allowing each molecular orbital to be a general linear combination of alpha and beta spin functions without any spin projection constraints. In contrast, Restricted Hartree–Fock (RHF) constrains electrons of alpha and beta spin to share identical spatial orbitals, which is suitable for closed-shell singlet systems but fails to describe bond dissociation. GHF relaxes these constraints, permitting the spin axis to vary spatially, which is essential for describing non-collinear magnetic structures and spin-orbit coupling. The GHF Roothaan–Hall equations are formulated in a combined basis of twice the dimension as

$$
\begin{pmatrix}\mathbf{F}^{\alpha\alpha}&\mathbf{F}^{\alpha\beta}\\\mathbf{F}^{\beta\alpha}&\mathbf{F}^{\beta\beta}\end{pmatrix}\begin{pmatrix}\mathbf{C}^{\alpha}\\\mathbf{C}^{\beta}\end{pmatrix}=\begin{pmatrix}\mathbf{S}&\mathbf{0}\\\mathbf{0}&\mathbf{S}\end{pmatrix}\begin{pmatrix}\mathbf{C}^{\alpha}\\\mathbf{C}^{\beta}\end{pmatrix}\mathbf{E}
$$

where the Fock matrix blocks $\mathbf{F}^{\alpha\alpha}$ and $\mathbf{F}^{\beta\beta}$ represent the spin-conserving Fock operators for alpha and beta spins, while the off-diagonal blocks $\mathbf{F}^{\alpha\beta}$ and $\mathbf{F}^{\beta\alpha}$ represent spin-mixing Fock operators that couple the alpha and beta spin components. The sub-matrices $\mathbf{C}^{\alpha}$ and $\mathbf{C}^{\beta}$ describe the spatial distribution of the alpha and beta spin components of the molecular orbitals, allowing spin polarization and non-collinear structures to be modeled.

---

## II. Energy and SCF Optimization

### 3. Fock Matrix Construction

For a restricted closed-shell system, the elements of the Fock matrix in the atomic orbital basis are constructed as

$$
F_{\mu\nu}=H_{\mu\nu}^{\text{core}}+\sum_{\lambda,\sigma}P_{\lambda\sigma}\left(\langle\mu\lambda|\nu\sigma\rangle-\frac{1}{2}\langle\mu\lambda|\sigma\nu\rangle\right)
$$

where $H_{\mu\nu}^{\text{core}}$ represents the core Hamiltonian matrix elements describing the kinetic energy of electrons and the electrostatic attraction between electrons and nuclei, $P_{\lambda\sigma}$ represents the density matrix elements, $\langle\mu\lambda|\nu\sigma\rangle$ is the two-electron Coulomb integral representing the classical electrostatic repulsion between electron densities, and $\langle\mu\lambda|\sigma\nu\rangle$ is the two-electron exchange integral representing the quantum mechanical correlation arising from the antisymmetrization of the wavefunction. The factor of $1/2$ on the exchange term is a spin-integration factor arising because the summation runs over spatial orbitals where each orbital is occupied by two electrons of opposite spin. The density matrix elements are computed from the occupied molecular orbital coefficients as

$$
P_{\lambda\sigma}=2\sum_i^{\text{occ}}C_{\lambda i}C_{\sigma i}
$$

where the factor of two accounts for the double occupancy of each spatial orbital $i$, and the summation runs over all occupied spatial molecular orbitals to define the total ground-state electronic density. In generalized Hartree–Fock, the Fock matrix and density matrix are expressed directly in terms of spin-orbitals, which avoids the spin-integration factors of two and yields a unified formulation where the exchange term directly cancels the self-interaction in the Coulomb term.

### 4. Self-Consistent Field Iteration and DIIS

Since the Fock matrix $\mathbf{F}$ depends on the density matrix $\mathbf{P}$, the Roothaan–Hall equations are nonlinear and must be solved iteratively. The self-consistent field cycle is accelerated using the Direct Inversion in the Iterative Subspace (DIIS) method, which minimizes the error matrix defined by the commutator of the Fock and density matrices in the overlap representation as

$$
\mathbf{e}=\mathbf{F}\mathbf{P}\mathbf{S}-\mathbf{S}\mathbf{P}\mathbf{F}
$$

where the error matrix $\mathbf{e}$ must vanish at convergence, reflecting that the Fock matrix and density matrix commute in the orthogonalized basis. The codebase implements the DIIS solver by building a subspace of historical Fock matrices and solving the linear system of the Lagrange multiplier constraints to obtain optimized interpolation weights.

### 5. Total Energy

The total electronic Hartree–Fock energy is calculated from the converged density and Fock matrices as

$$
E_{\text{elec}}=\frac{1}{2}\sum_{\mu,\nu}P_{\mu\nu}\left(H_{\mu\nu}^{\text{core}}+F_{\mu\nu}\right)
$$

where the factor of $1/2$ prevents the double-counting of electron-electron repulsion energy contained within the Fock matrix elements, and the total energy of the system is the sum of the electronic energy and the classical nuclear repulsion energy.

---

## III. Nuclear Derivatives and Response

### 6. Analytical Nuclear Gradient

The analytical gradient of the Hartree–Fock energy with respect to a nuclear coordinate $x$ is calculated using the density matrix and the derivatives of the integrals as

$$
\frac{dE_{\text{HF}}}{dx}=\frac{dV_{\text{nuc}}}{dx}+\sum_{\mu,\nu}P_{\mu\nu}\frac{dH_{\mu\nu}^{\text{core}}}{dx}-\sum_{\mu,\nu}W_{\mu\nu}\frac{dS_{\mu\nu}}{dx}+\frac{1}{2}\sum_{\mu,\nu,\lambda,\sigma}P_{\mu\nu}P_{\lambda\sigma}\left(\frac{d\langle\mu\lambda|\nu\sigma\rangle}{dx}-\frac{1}{2}\frac{d\langle\mu\lambda|\sigma\nu\rangle}{dx}\right)
$$

where $V_{\text{nuc}}$ is the classical nuclear repulsion energy, $H_{\mu\nu}^{\text{core}}$ is the core Hamiltonian, $S_{\mu\nu}$ is the overlap matrix, $P_{\mu\nu}$ is the density matrix, and $W_{\mu\nu}$ is the energy-weighted density matrix defined as

$$
W_{\mu\nu}=2\sum_i^{\text{occ}}\epsilon_iC_{\mu i}C_{\nu i}
$$

with $\epsilon_i$ representing the molecular orbital energies. The derivatives of the molecular orbital coefficients do not appear explicitly in the gradient expression due to the variational optimization of the Hartree–Fock wavefunction, meaning the energy is stationary with respect to orbital variations; the overlap derivative term represents the Pulay force arising because the atomic orbital basis functions are atom-centered and move with the nuclei.

### 7. Coupled-Perturbed Hartree–Fock Equations

The Coupled-Perturbed Hartree–Fock (CPHF) equations describe the linear response of the molecular orbital coefficients to a perturbation such as a nuclear displacement or an external electric field. Because the Fock matrix depends on the density matrix, any change in the molecular orbitals feeds back into the Fock matrix, requiring a self-consistent determination of the orbital response. We express the derivative of the molecular orbital coefficients in terms of the unperturbed coefficients using the orbital response matrix $\mathbf{U}^x$ as

$$
\frac{dC_{\mu i}}{dx}=\sum_pC_{\mu p}U_{pi}^x
$$

where $p$ runs over all occupied and virtual molecular orbitals, defining how the unperturbed molecular orbitals mix under the perturbation. The orthonormality constraint of the molecular orbitals determines the symmetric part of the response matrix, which relates to the derivative of the overlap matrix $\mathbf{S}$ in the molecular orbital basis as

$$
U_{pq}^x+U_{qp}^x+\frac{dS_{pq}}{dx}=0
$$

where $S_{pq}$ is the overlap matrix in the molecular orbital basis, ensuring the perturbed orbitals remain orthonormal. The remaining occupied-virtual blocks of the response matrix are obtained by solving the linear CPHF equations

$$
(\epsilon_a-\epsilon_i)U_{ai}^x-\sum_{j}^{\text{occ}}\sum_{b}^{\text{vir}}A_{ai,bj}U_{bj}^x=B_{ai}^x
$$

where $i$ and $j$ denote occupied orbitals, $a$ and $b$ denote virtual orbitals, $\epsilon_p$ represents the orbital energies, $A_{ai,bj}$ are the coupling matrix elements reflecting the change in the Hartree–Fock potential due to density matrix modifications, and $B_{ai}^x$ represents the direct perturbation vector containing the derivatives of the Fock and overlap matrices in the molecular orbital basis. To avoid explicitly constructing and inverting the large coupling matrix $\mathbf{A}$, the codebase implements an iterative CPHF solver that converges the response equations using DIIS subspace acceleration.
