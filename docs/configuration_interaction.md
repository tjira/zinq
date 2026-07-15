# Configuration–Interaction

The Configuration–Interaction (CI) method is a linear variational post-Hartree–Fock method to describe electron correlation in many-electron systems. It expands the N-electron wavefunction as a linear combination of Slater determinants corresponding to different electronic configurations, including the reference Hartree–Fock determinant and various excited configurations.

---

## I. Wavefunction Representation

### 1. Variational Wavefunction Expansion

The CI wavefunction is expanded as

$$
| \Psi_{\text{CI}} \rangle = c_0 | \Phi_0 \rangle + \sum_i c_i | \Phi_i \rangle
$$

where $| \Phi_0 \rangle$ is the reference Hartree–Fock Slater determinant representing the mean-field ground state, $| \Phi_i \rangle$ represents the excited Slater determinants where one or more electrons are promoted from occupied to virtual spatial/spin orbitals, and $c_i$ are the variational coefficients. In determinant-based Configuration–Interaction, the wavefunction expansion uses Slater determinants directly, which simplifies the evaluation of matrix elements but does not guarantee that the resulting wavefunction is an eigenstate of the total spin operator $\hat{S}^2$. The codebase represents Slater determinants using arrays of occupied spin-orbital indices, permitting direct matching and excitation degree evaluation.

---

## II. Hamiltonian Evaluation and Diagonalization

### 2. Slater–Condon Rules

The matrix elements $H_{ij}=\langle\Phi_i|\hat{H}|\Phi_j\rangle$ of the electronic Hamiltonian between Slater determinants are evaluated using the Slater–Condon rules. These rules simplify the matrix elements based on the degree of excitation between the two determinants, reducing many-electron integrals to one- and two-electron molecular integrals under the assumption that the spin-orbitals are orthonormal. If the two determinants are identical, the matrix elements of a one-body operator $\hat{F}=\sum_k\hat{f}(k)$ and a two-body operator $\hat{G}=\sum_{k<l}\hat{g}(k,l)$ are given by

$$
\langle\Phi_A|\hat{F}|\Phi_A\rangle=\sum_i\langle\phi_i|\hat{f}|\phi_i\rangle
$$

and

$$
\langle\Phi_A|\hat{G}|\Phi_A\rangle=\frac{1}{2}\sum_{i,j}\langle\phi_i\phi_j||\phi_i\phi_j\rangle
$$

where $\hat{f}(k)$ and $\hat{g}(k,l)$ are the one- and two-body operators acting on electrons $k$ and $l$ respectively, the sums run over all occupied spin-orbitals $\phi_i$, and the double bar indicates an antisymmetrized two-electron integral $\langle ij || ij \rangle = \langle ij | ij \rangle - \langle ij | ji \rangle$ combining the Coulomb and exchange interactions. For determinants differing by a single spin-orbital substitution from $\phi_p$ to $\phi_r$, the matrix elements are

$$
\langle\Phi_A|\hat{F}|\Phi_B\rangle=\langle\phi_p|\hat{f}|\phi_r\rangle
$$

and

$$
\langle\Phi_A|\hat{G}|\Phi_B\rangle=\sum_i\langle\phi_p\phi_i||\phi_r\phi_i\rangle
$$

where the sum runs over the occupied spin-orbitals common to both determinants, representing the transition density interaction. For determinants differing by two spin-orbital substitutions from $\phi_p,\phi_q$ to $\phi_r,\phi_s$, the one-body matrix element vanishes and the two-body matrix element is

$$
\langle\Phi_A|\hat{G}|\Phi_B\rangle=\langle\phi_p\phi_q||\phi_r\phi_s\rangle
$$

while all matrix elements vanish if the determinants differ by three or more spin-orbitals because the Hamiltonian only contains up to two-body interactions. To evaluate these elements, the codebase compares the occupied spin-orbital indices of the two determinants to determine the excitation degree. The phase sign of the excitation is computed by counting the number of transpositions (permutations) needed to align the spin-orbital indices, ensuring the antisymmetry of the electronic wavefunction is preserved exactly.

### 3. Hamiltonian Diagonalization

Applying the variational principle to optimize the coefficients leads to the matrix eigenvalue equation

$$
\mathbf{H}_{\text{CI}}\mathbf{C}_k=E_k\mathbf{C}_k
$$

where $\mathbf{H}_{\text{CI}}$ is the Hamiltonian matrix representation in the Slater determinant basis, $\mathbf{C}_k$ is the eigenvector containing the coefficients for state $k$, and $E_k$ is the corresponding electronic energy. The total energy of state $k$ is the sum of the electronic energy and the classical nuclear repulsion energy. The codebase solves this generalized eigenvalue problem using symmetric matrix diagonalization from `linear_algebra.zig`.

---

## III. Nuclear Derivatives

### 4. Analytical Nuclear Gradient

The analytical gradient of the Configuration–Interaction energy of state $k$ with respect to a nuclear coordinate $x$ is formulated using the Hellmann–Feynman theorem. Since the CI wavefunction coefficients are variationally optimized by diagonalizing the Hamiltonian matrix, their derivatives with respect to the nuclear coordinates do not contribute to the energy derivative. The resulting analytical nuclear gradient is given by

$$
\frac{dE_k}{dx}=\frac{dV_{\text{nuc}}}{dx}+\sum_{a,b}C_{ak}C_{bk}\frac{dH_{\text{CI},ab}}{dx}
$$

where $V_{\text{nuc}}$ is the classical nuclear repulsion energy, $C_{ak}$ and $C_{bk}$ are the components of the CI eigenvector $\mathbf{C}_k$ for state $k$ on determinants $a$ and $b$, and $H_{\text{CI},ab}$ is the CI Hamiltonian matrix $\mathbf{H}_{\text{CI}}$ elements in the Slater determinant basis. The derivatives of the Hamiltonian matrix elements $dH_{\text{CI},ab}/dx$ are evaluated by propagating the nuclear derivatives of the core Hamiltonian, overlap, and two-electron integrals, alongside the molecular orbital coefficient derivatives obtained from coupled-perturbed Hartree–Fock (CPHF) equations to account for the implicit coordinate dependence of the molecular orbitals.
