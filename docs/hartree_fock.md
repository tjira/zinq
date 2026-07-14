# Hartree–Fock Method

The Hartree–Fock method is a fundamental variational approach in quantum chemistry and electronic structure theory to approximate the ground-state electronic wavefunction of a multi-electron system. It models the wavefunction as a single Slater determinant of molecular orbitals, effectively introducing a mean-field approximation where each electron moves in the average electrostatic field of all other electrons.

---

## 1. Roothaan–Hall Equations

By expanding the molecular orbitals as a linear combination of atomic orbitals, the variational optimization of the energy leads to the Roothaan–Hall equations

$$
\mathbf{F}\mathbf{C}=\mathbf{S}\mathbf{C}\mathbf{E}
$$

where $\mathbf{F}$ is the Fock matrix, $\mathbf{C}$ is the molecular orbital coefficient matrix, $\mathbf{S}$ is the atomic orbital overlap matrix, and $\mathbf{E}$ is the diagonal matrix of molecular orbital energies.

---

## 2. Fock Matrix Construction

For a restricted closed-shell system, the elements of the Fock matrix in the atomic orbital basis are constructed as

$$
F_{\mu\nu}=H_{\mu\nu}^{\text{core}}+\sum_{\lambda,\sigma}P_{\lambda\sigma}\left((\mu\nu||\lambda\sigma)-\frac{1}{2}(\mu\lambda||\nu\sigma)\right)
$$

where $H_{\mu\nu}^{\text{core}}$ represents the core Hamiltonian matrix elements describing kinetic energy and electron-nuclear attraction, $P_{\lambda\sigma}$ represents the density matrix elements, $(\mu\nu||\lambda\sigma)$ is the two-electron Coulomb integral, and $(\mu\lambda||\nu\sigma)$ is the two-electron exchange integral. The density matrix elements are computed from the occupied molecular orbital coefficients as

$$
P_{\lambda\sigma}=2\sum_i^{\text{occ}}C_{\lambda i}C_{\sigma i}
$$

which defines the electronic density of the system.

---

## 3. Self-Consistent Field Iteration and DIIS

Since the Fock matrix $\mathbf{F}$ depends on the density matrix $\mathbf{P}$, the Roothaan–Hall equations are nonlinear and must be solved iteratively. The self-consistent field cycle is accelerated using the Direct Inversion in the Iterative Subspace (DIIS) method, which minimizes the error matrix defined by the commutator of the Fock and density matrices in the overlap representation as

$$
\mathbf{e}=\mathbf{F}\mathbf{P}\mathbf{S}-\mathbf{S}\mathbf{P}\mathbf{F}
$$

where the convergence is monitored by the magnitude of the error commutator.

---

## 4. Total Energy

The total electronic Hartree–Fock energy is calculated from the converged density and Fock matrices as

$$
E_{\text{elec}}=\frac{1}{2}\sum_{\mu,\nu}P_{\mu\nu}\left(H_{\mu\nu}^{\text{core}}+F_{\mu\nu}\right)
$$

and the total energy of the system is the sum of the electronic energy and the classical nuclear repulsion energy.

---

## 5. Analytical Nuclear Gradient

The analytical gradient of the Hartree–Fock energy with respect to a nuclear coordinate $x$ is calculated using the density matrix and the derivatives of the integrals as

$$
\frac{dE_{\text{HF}}}{dx}=\frac{dV_{\text{nuc}}}{dx}+\sum_{\mu,\nu}P_{\mu\nu}\frac{dH_{\mu\nu}^{\text{core}}}{dx}-\sum_{\mu,\nu}W_{\mu\nu}\frac{dS_{\mu\nu}}{dx}+\frac{1}{2}\sum_{\mu,\nu,\lambda,\sigma}P_{\mu\nu}P_{\lambda\sigma}\left(\frac{d(\mu\nu||\lambda\sigma)}{dx}-\frac{1}{2}\frac{d(\mu\lambda||\nu\sigma)}{dx}\right)
$$

where $V_{\text{nuc}}$ is the classical nuclear repulsion energy, $H_{\mu\nu}^{\text{core}}$ is the core Hamiltonian, $S_{\mu\nu}$ is the overlap matrix, $P_{\mu\nu}$ is the density matrix, and $W_{\mu\nu}$ is the energy-weighted density matrix defined as

$$
W_{\mu\nu}=2\sum_i^{\text{occ}}\epsilon_iC_{\mu i}C_{\nu i}
$$

with $\epsilon_i$ representing the orbital energies. The derivatives of the molecular orbital coefficients do not appear explicitly in the gradient expression due to the variational optimization of the Hartree–Fock wavefunction.
