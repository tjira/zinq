# Configuration–Interaction

The Configuration–Interaction (CI) method is a linear variational post-Hartree–Fock method to describe electron correlation in many-electron systems. It expands the N-electron wavefunction as a linear combination of Slater determinants corresponding to different electronic configurations, including the reference Hartree–Fock determinant and various excited configurations.

---

## 1. Variational Wavefunction Expansion

The CI wavefunction is expanded as

$$
| \Psi_{\text{CI}} \rangle = c_0 | \Phi_0 \rangle + \sum_i c_i | \Phi_i \rangle
$$

where $| \Phi_0 \rangle$ is the reference Hartree–Fock Slater determinant, $| \Phi_i \rangle$ represents the excited Slater determinants, and $c_i$ are the variational coefficients.

---

## 2. Hamiltonian Diagonalization

Applying the variational principle to optimize the coefficients leads to the matrix eigenvalue equation

$$
\mathbf{H}_{\text{CI}}\mathbf{C}_k=E_k\mathbf{C}_k
$$

where $\mathbf{H}_{\text{CI}}$ is the Hamiltonian matrix representation in the Slater determinant basis, $\mathbf{C}_k$ is the eigenvector containing the coefficients for state $k$, and $E_k$ is the corresponding electronic energy. The total energy of state $k$ is the sum of the electronic energy and the classical nuclear repulsion energy.

---

## 3. Slater–Condon Rules

The matrix elements $H_{ij}=\langle\Phi_i|\hat{H}|\Phi_j\rangle$ of the electronic Hamiltonian between Slater determinants are evaluated using the Slater–Condon rules. These rules simplify the matrix elements based on the degree of excitation between the two determinants, reducing the N-electron integrals to one- and two-electron molecular integrals.

---

## 4. Analytical Nuclear Gradient

The analytical gradient of the Configuration–Interaction energy of state $k$ with respect to a nuclear coordinate $x$ is formulated using the Hellmann–Feynman theorem. Since the CI wavefunction coefficients are variationally optimized by diagonalizing the Hamiltonian matrix, their derivatives with respect to the nuclear coordinates do not contribute to the energy derivative. The resulting analytical nuclear gradient is given by

$$
\frac{dE_k}{dx}=\frac{dV_{\text{nuc}}}{dx}+\sum_{a,b}C_{ak}C_{bk}\frac{dH_{\text{CI},ab}}{dx}
$$

where $V_{\text{nuc}}$ is the classical nuclear repulsion energy, $C_{ak}$ and $C_{bk}$ are the components of the CI eigenvector $\mathbf{C}_k$ for state $k$ on determinants $a$ and $b$, and $H_{\text{CI},ab}$ is the CI Hamiltonian matrix $\mathbf{H}_{\text{CI}}$ elements in the Slater determinant basis. The derivatives of the Hamiltonian matrix elements $dH_{\text{CI},ab}/dx$ are evaluated by propagating the nuclear derivatives of the core Hamiltonian, overlap, and two-electron integrals, alongside the molecular orbital coefficient derivatives obtained from coupled-perturbed Hartree–Fock (CPHF) equations.
