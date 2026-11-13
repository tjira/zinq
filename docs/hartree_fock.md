# Hartree–Fock Method

This document provides a clear, simple, and mathematically rigorous explanation of the Hartree–Fock method in second quantization. If you want to understand how the electronic Hamiltonian, Slater determinants, Wick's theorem, basis sets, the variational principle, and Brillouin's theorem lead to the fundamental mean-field equations of quantum chemistry, this guide explains every step directly.

---

## I. Second-Quantized Hamiltonian and Reference State

### 1. Exact Electronic Hamiltonian in Second Quantization

In second quantization, we describe electrons using creation and annihilation operators. This makes it easy to keep track of the antisymmetry of electrons without having to write out large determinants by hand. We begin by choosing an orthonormal set of spin-orbitals $\{\phi_p(\mathbf{x})\}$, where $\mathbf{x}=(\mathbf{r},\sigma)$ includes both spatial coordinates $\mathbf{r}$ and spin coordinates $\sigma$. The creation operator $a_p^\dagger$ adds an electron to orbital $p$, and the annihilation operator $a_q$ removes an electron from orbital $q$. Because electrons are fermions, these operators satisfy canonical anticommutation relations in square brackets as

$$
[a_p,a_q^\dagger]_+=\delta_{pq},\quad[a_p,a_q]_+=[a_p^\dagger,a_q^\dagger]_+=0
$$

where $[A,B]_+=AB+BA$, and $\delta_{pq}$ is the Kronecker delta. The exact non-relativistic electronic Hamiltonian operator $\hat{H}$ is written with one-electron and two-electron terms using these operators as

$$
\hat{H}=\sum_{pq}h_{pq}a_p^\dagger a_q+\frac{1}{2}\sum_{pqrs}g_{pqrs}a_p^\dagger a_q^\dagger a_s a_r
$$

where $h_{pq}=\langle\phi_p|\hat{h}|\phi_q\rangle$ represents the one-electron core Hamiltonian describing the kinetic energy of an electron and its electrostatic attraction to all nuclei of charge $Z_A$ at positions $\mathbf{R}_A$, with the one-electron operator $\hat{h}$ defined as

$$
\hat{h}=-\frac{1}{2}\nabla^2-\sum_A\frac{Z_A}{|\mathbf{r}-\mathbf{R}_A|}
$$

and $g_{pqrs}=\langle pq|rs\rangle$ represents the two-electron repulsion integral in physicist notation describing the electrostatic repulsion between two electrons as

$$
g_{pqrs}=\langle pq|rs\rangle=\int\int\phi_p^*(\mathbf{x}_1)\phi_q^*(\mathbf{x}_2)\frac{1}{r_{12}}\phi_r(\mathbf{x}_1)\phi_s(\mathbf{x}_2)d\mathbf{x}_1d\mathbf{x}_2
$$

with $r_{12}=|\mathbf{r}_1-\mathbf{r}_2|$. The factor of $1/2$ avoids double-counting each electron pair repulsion.

### 2. Hartree–Fock Reference State and Wick's Theorem Evaluation

Next, we define the Hartree–Fock reference state $|\Phi_0\rangle$ as a single Slater determinant of $N$ occupied electrons. In second quantization, we create this state by applying creation operators for each occupied orbital to the empty vacuum state $|0\rangle$ as

$$
|\Phi_0\rangle=\prod_{i=1}^N a_i^\dagger|0\rangle=a_1^\dagger a_2^\dagger\cdots a_N^\dagger|0\rangle
$$

where indices $i,j,k,\dots$ label occupied orbitals and $a,b,c,\dots$ label empty virtual orbitals. The state $|\Phi_0\rangle$ serves as the Fermi vacuum of the system. In this picture, occupied orbitals act as holes and empty orbitals act as particles, meaning that trying to annihilate an electron in an unoccupied virtual orbital or create an electron in an already occupied orbital vanishes by the Pauli exclusion principle, so $a_a|\Phi_0\rangle=0$ and $a_i^\dagger|\Phi_0\rangle=0$.

To find the energy of this state, we evaluate the expectation value of the Hamiltonian with respect to the Fermi vacuum, given by

$$
E_{\text{HF}}=\langle\Phi_0|\hat{H}|\Phi_0\rangle=\sum_{pq}h_{pq}\langle\Phi_0|a_p^\dagger a_q|\Phi_0\rangle+\frac{1}{2}\sum_{pqrs}g_{pqrs}\langle\Phi_0|a_p^\dagger a_q^\dagger a_s a_r|\Phi_0\rangle
$$

To calculate the expectation values of these operator strings, we use Wick's theorem for expectation values. For the Fermi vacuum $|\Phi_0\rangle$, Wick's theorem states that the expectation value of any string of creation and annihilation operators equals the sum of all possible fully contracted pairs as

$$
\langle\Phi_0|A_1 A_2\cdots A_{2m}|\Phi_0\rangle=\sum_{\text{pairings}}(-1)^\sigma\prod_{k=1}^m\langle\Phi_0|A_{i_k}A_{j_k}|\Phi_0\rangle
$$

where the sum runs over all ways to partition the $2m$ operators into $m$ disjoint pairs $(A_{i_k},A_{j_k})$ with position indices satisfying $i_k<j_k$, so that $i_k$ and $j_k$ are the positions of the two operators contracted together in pair $k$, and $(-1)^\sigma$ is the sign of the permutation that rearranges the original sequence $(1,2,\dots,2m)$ into the paired sequence $(i_1,j_1,i_2,j_2,\dots,i_m,j_m)$. If a string has an odd number of operators, its expectation value is zero.

Because the Hamiltonian is normal-ordered with respect to the bare vacuum, every contracted pair places a creation operator to the left of an annihilation operator. For the Fermi vacuum reference state $|\Phi_0\rangle$, annihilation operators kill virtual orbitals to the right, $a_a|\Phi_0\rangle=0$, while creation operators kill virtual orbitals to the left on the dual bra, $\langle\Phi_0|a_a^\dagger=0$. The only contracted pair with a non-zero expectation value therefore occurs when both operators belong to occupied spin-orbitals as

$$
\langle\Phi_0|a_p^\dagger a_q|\Phi_0\rangle=\delta_{pq}\quad(p,q\in\text{occ})
$$

while any pair containing an unoccupied virtual orbital vanishes ($p\in\text{virt}$ or $q\in\text{virt}$). Because each spin-orbital in a single determinant is either strictly occupied or strictly empty, occupation numbers are not needed. Applying Wick's theorem to the one-body term therefore eliminates all virtual orbitals and simplifies the one-body energy to $\sum_i^{\text{occ}}h_{ii}$.

Similarly, for the two-body operator string $a_p^\dagger a_q^\dagger a_s a_r$, the expectation value is non-zero only when all four indices belong to occupied spin-orbitals ($p,q,r,s\in\text{occ}$). If any annihilation operator belongs to a virtual orbital ($r\text{ or }s\in\text{virt}$), it destroys the ket $|\Phi_0\rangle$, while if any creation operator belongs to a virtual orbital ($p\text{ or }q\in\text{virt}$), its adjoint destroys the bra $\langle\Phi_0|$. Applying Wick's theorem to the string of occupied spin-orbitals produces two non-zero pairings as

$$
\langle\Phi_0|a_p^\dagger a_q^\dagger a_s a_r|\Phi_0\rangle=\langle\Phi_0|a_p^\dagger a_r|\Phi_0\rangle\langle\Phi_0|a_q^\dagger a_s|\Phi_0\rangle-\langle\Phi_0|a_p^\dagger a_s|\Phi_0\rangle\langle\Phi_0|a_q^\dagger a_r|\Phi_0\rangle=\delta_{pr}\delta_{qs}-\delta_{ps}\delta_{qr}\quad(p,q,r,s\in\text{occ})
$$

where the minus sign comes from swapping two fermionic operators. Contracting the creation and annihilation operators eliminates all terms except those involving occupied spin-orbitals, reducing the two-body sum to $\frac{1}{2}\sum_{i,j}^{\text{occ}}\left(g_{ijij}-g_{ijji}\right)$, where $g_{ijij}$ is the direct Coulomb repulsion and $g_{ijji}$ is the quantum exchange interaction. The total Hartree–Fock energy is therefore given by the algebraic expression

$$
E_{\text{HF}}=\sum_i^{\text{occ}}h_{ii}+\frac{1}{2}\sum_{i,j}^{\text{occ}}\left(g_{ijij}-g_{ijji}\right)
$$

which depends only on the core Hamiltonian elements $h_{ii}$ and the two-electron integrals $g_{pqrs}$ for the occupied orbitals.

---

## II. Algebraic Approximation and Basis Set Expansion

### 3. Roothaan–Hall Expansion and the One-Particle Density Matrix

Although this energy expression is exact for a single determinant, solving the spatial differential equations directly in three dimensions is intractable for general molecules because the potential from multiple nuclei is not spherical. To solve these equations in practice, we introduce the algebraic approximation by using a known basis set. In Dirac notation, Generalized Hartree–Fock (GHF) expands each unknown molecular spin-orbital state $|\phi_p\rangle$ as a linear combination of known, non-orthogonal spin-orbital basis states $|\chi_\mu\rangle$ centered on the atoms as

$$
|\phi_p\rangle=\sum_\mu C_{\mu p}|\chi_\mu\rangle
$$

where $C_{\mu p}$ are the complex expansion coefficients we want to find ($C_{\mu p}\in\mathbb{C}$), and the corresponding dual bra is $\langle\phi_p|=\sum_\mu C_{\mu p}^*\langle\chi_\mu|$. The requirement that molecular orbitals stay orthonormal, $\langle\phi_p|\phi_q\rangle=\delta_{pq}$, is evaluated by inserting the expansion for both the bra and the ket directly into the Dirac bracket as

$$
\langle\phi_p|\phi_q\rangle=\left(\sum_\mu C_{\mu p}^*\langle\chi_\mu|\right)\left(\sum_\nu C_{\nu q}|\chi_\nu\rangle\right)=\sum_{\mu,\nu}C_{\mu p}^*\langle\chi_\mu|\chi_\nu\rangle C_{\nu q}=\sum_{\mu,\nu}C_{\mu p}^*S_{\mu\nu}C_{\nu q}=\delta_{pq}
$$

where $S_{\mu\nu}=\langle\chi_\mu|\chi_\nu\rangle$ is the atomic orbital overlap matrix element measuring the spatial overlap between the non-orthogonal basis functions. In matrix form, this orthonormality condition is written compactly as $\mathbf{C}^\dagger\mathbf{S}\mathbf{C}=\mathbf{I}$, where $\mathbf{I}$ is the identity matrix and $\mathbf{C}^\dagger$ is the conjugate transpose of $\mathbf{C}$. Because atomic orbitals centered on different nuclei are not mutually orthogonal ($S_{\mu\nu}\neq\delta_{\mu\nu}$), the overlap matrix $\mathbf{S}$ acts as a metric ensuring that the molecular orbitals remain orthonormal in space.

To represent the electrons in this atomic basis, we introduce the one-particle density matrix $\mathbf{P}$ by summing over all occupied molecular orbitals as

$$
P_{\mu\nu}=\sum_i^{\text{occ}}C_{\mu i}C_{\nu i}^*
$$

which directly stores the spatial distribution of the electrons in the atomic orbital basis, satisfying the Hermitian symmetry relation $P_{\mu\nu}=P_{\nu\mu}^*$.

### 4. Projection onto the Atomic Orbital Basis

Next, we project the energy expression onto the atomic orbital basis using this Dirac representation. For the one-electron core energy, the expectation value of the one-body operator $\hat{h}$ for an occupied orbital $|\phi_i\rangle$ is projected by inserting the basis expansions as

$$
h_{ii}=\langle\phi_i|\hat{h}|\phi_i\rangle=\left(\sum_\mu C_{\mu i}^*\langle\chi_\mu|\right)\hat{h}\left(\sum_\nu C_{\nu i}|\chi_\nu\rangle\right)=\sum_{\mu,\nu}C_{\mu i}^*\langle\chi_\mu|\hat{h}|\chi_\nu\rangle C_{\nu i}
$$

Because $\hat{h}$ is linear, the coefficients pull outside the bracket, leaving the atomic core Hamiltonian matrix element $H_{\mu\nu}=\langle\chi_\mu|\hat{h}|\chi_\nu\rangle$. This matrix element depends only on the predefined atomic basis states and the nuclear potential, so it can be computed once and stored as a fixed matrix $\mathbf{H}$. Summing over all occupied molecular orbitals $i$ groups the coefficients into the Hermitian density matrix $P_{\nu\mu}=\sum_i^{\text{occ}}C_{\nu i}C_{\mu i}^*$, giving the total one-electron core energy as

$$
\sum_i^{\text{occ}}h_{ii}=\sum_{\mu,\nu}P_{\nu\mu}H_{\mu\nu}
$$

We project the two-electron repulsion integrals $g_{ijij}-g_{ijji}$ onto the atomic basis in the exact same way. In physicist notation, the atomic two-electron integral is written as the Dirac bracket

$$
\langle\mu\lambda|\nu\sigma\rangle=\int\int\chi_\mu^*(\mathbf{x}_1)\chi_\lambda^*(\mathbf{x}_2)\frac{1}{r_{12}}\chi_\nu(\mathbf{x}_1)\chi_\sigma(\mathbf{x}_2)d\mathbf{x}_1d\mathbf{x}_2
$$

where electron 1 occupies basis states $\mu$ and $\nu$, while electron 2 occupies basis states $\lambda$ and $\sigma$. To see how the contraction with the density matrix arises, we expand the molecular orbitals in both the Coulomb and exchange terms. For the direct Coulomb term $g_{ijij}=\langle ij|ij\rangle$, inserting the expansions for all four orbitals yields

$$
\langle ij|ij\rangle=\sum_{\mu,\nu,\lambda,\sigma}C_{\mu i}^*C_{\nu i}C_{\lambda j}^*C_{\sigma j}\langle\mu\lambda|\nu\sigma\rangle
$$

Summing over all occupied orbitals $i$ and $j$ separates the sums into independent products of coefficients. Because $\sum_i^{\text{occ}}C_{\mu i}^*C_{\nu i}$ and $\sum_j^{\text{occ}}C_{\lambda j}^*C_{\sigma j}$ contract with the symmetric two-electron integrals into density matrices, the total Coulomb repulsion simplifies to

$$
\sum_{i,j}^{\text{occ}}g_{ijij}=\sum_{\mu,\nu,\lambda,\sigma}P_{\mu\nu}P_{\lambda\sigma}\langle\mu\lambda|\nu\sigma\rangle
$$

For the exchange term $g_{ijji}=\langle ij|ji\rangle$, the ket orbitals are swapped, meaning electron 1 is in orbital $\phi_j$ and electron 2 is in orbital $\phi_i$. Inserting the basis expansions swaps the ket basis states to $\langle\mu\lambda|\sigma\nu\rangle$ as

$$
\langle ij|ji\rangle=\sum_{\mu,\nu,\lambda,\sigma}C_{\mu i}^*C_{\nu i}C_{\lambda j}^*C_{\sigma j}\langle\mu\lambda|\sigma\nu\rangle
$$

Summing over occupied orbitals $i$ and $j$ contracts these coefficients into the exact same density matrices $P_{\mu\nu}$ and $P_{\lambda\sigma}$, giving the total exchange energy as

$$
\sum_{i,j}^{\text{occ}}g_{ijji}=\sum_{\mu,\nu,\lambda,\sigma}P_{\mu\nu}P_{\lambda\sigma}\langle\mu\lambda|\sigma\nu\rangle
$$

Subtracting the exchange contribution from the Coulomb contribution with the prefactor $1/2$ and adding the one-electron core energy yields the total Hartree–Fock energy in closed form as

$$
E_{\text{HF}}=\sum_{\mu,\nu}P_{\nu\mu}H_{\mu\nu}+\frac{1}{2}\sum_{\mu,\nu,\lambda,\sigma}P_{\mu\nu}P_{\lambda\sigma}\left(\langle\mu\lambda|\nu\sigma\rangle-\langle\mu\lambda|\sigma\nu\rangle\right)
$$

which depends only on the core Hamiltonian matrix $\mathbf{H}$, the density matrix $\mathbf{P}$, and the two-electron repulsion integrals $\langle\mu\lambda|\nu\sigma\rangle$.

---

## III. Variational Optimization and the Roothaan–Hall Equations

### 5. Variational Principle and Fock Operator Derivation

With the energy expressed in a finite basis, we apply the variational principle to find the electronic ground state. The variational principle states that the optimal orbital coefficients are the ones that minimize the total energy $E_{\text{HF}}$, while keeping the molecular orbitals orthonormal, $\sum_{\mu,\nu}C_{\mu i}^*S_{\mu\nu}C_{\nu j}=\delta_{ij}$. To perform this constrained minimization, we construct the Lagrangian functional $\mathcal{L}$ with Lagrange multipliers $\epsilon_{ji}$ as

$$
\mathcal{L}=E_{\text{HF}}-\sum_{k,l}^{\text{occ}}\epsilon_{lk}\left(\sum_{\alpha,\beta}C_{\alpha k}^*S_{\alpha\beta}C_{\beta l}-\delta_{kl}\right)
$$

where the multipliers $\epsilon_{lk}$ enforce the orthonormality constraint for each pair of occupied orbitals $k$ and $l$.

In complex variation, we treat the dual bra coefficients $C_{\mu i}^*$ and ket coefficients $C_{\nu i}$ as independent variables using Wirtinger calculus, setting the variation with respect to the bra coefficients to zero, $\frac{\partial\mathcal{L}}{\partial C_{\mu i}^*}=0$. Taking the derivative of the constraint term gives

$$
\frac{\partial}{\partial C_{\mu i}^*}\left[-\sum_{k,l}^{\text{occ}}\epsilon_{lk}\sum_{\alpha,\beta}C_{\alpha k}^*S_{\alpha\beta}C_{\beta l}\right]=-\sum_{k,l}^{\text{occ}}\epsilon_{lk}\sum_{\alpha,\beta}\delta_{\alpha\mu}\delta_{ki}S_{\alpha\beta}C_{\beta l}=-\sum_j^{\text{occ}}\sum_\nu S_{\mu\nu}C_{\nu j}\epsilon_{ji}
$$

Next, we differentiate the energy $E_{\text{HF}}$ with respect to $C_{\mu i}^*$. The one-electron core term $\sum_{k}^{\text{occ}}\sum_{\alpha,\beta}C_{\alpha k}^*H_{\alpha\beta}C_{\beta k}$ gives $\sum_\nu H_{\mu\nu}C_{\nu i}$. Differentiating the two-electron term with respect to $C_{\mu i}^*$ acts on the bra of occupied orbital $i$, which appears symmetrically in both electron positions, cancelling the prefactor $1/2$ and yielding

$$
\frac{\partial E_{\text{HF}}}{\partial C_{\mu i}^*}=\sum_\nu H_{\mu\nu}C_{\nu i}+\sum_{\nu,\lambda,\sigma}P_{\lambda\sigma}\left(\langle\mu\lambda|\nu\sigma\rangle-\langle\mu\lambda|\sigma\nu\rangle\right)C_{\nu i}
$$

Combining the energy and constraint derivatives and setting $\frac{\partial\mathcal{L}}{\partial C_{\mu i}^*}=0$ gives the stationary condition

$$
\sum_\nu\left[H_{\mu\nu}+\sum_{\lambda,\sigma}P_{\lambda\sigma}\left(\langle\mu\lambda|\nu\sigma\rangle-\langle\mu\lambda|\sigma\nu\rangle\right)\right]C_{\nu i}=\sum_{j}^{\text{occ}}\sum_\nu S_{\mu\nu}C_{\nu j}\epsilon_{ji}
$$

The term in square brackets defines the matrix element $F_{\mu\nu}=\langle\chi_\mu|\hat{F}|\chi_\nu\rangle$ of an effective one-electron operator called the Fock operator $\hat{F}$. Recognizing that $\sum_\nu F_{\mu\nu}C_{\nu i}=\langle\chi_\mu|\hat{F}|\phi_i\rangle$ and $\sum_\nu S_{\mu\nu}C_{\nu j}=\langle\chi_\mu|\phi_j\rangle$, the stationary condition is the atomic basis representation of the operator equation

$$
\hat{F}|\phi_i\rangle=\sum_j^{\text{occ}}\epsilon_{ji}|\phi_j\rangle
$$

Minimizing the electronic energy under the orthonormality constraint replaces the interacting many-electron problem with an effective one-electron eigenvalue problem, where each electron moves in the field of the nuclei and the average mean-field potential generated by all electrons combined.

### 6. Coulomb, Exchange, and Roothaan–Hall Equations

The Fock operator $\hat{F}$ is the effective one-electron Hamiltonian describing an electron moving in the field of the nuclei and the mean field of all other electrons, expressed as

$$
\hat{F}=\hat{h}+\hat{V}^{\text{HF}}=\hat{h}+\hat{J}-\hat{K}
$$

where $\hat{h}$ is the one-electron core Hamiltonian representing kinetic energy and nuclear attraction, $\hat{J}$ is the classical Coulomb operator, and $\hat{K}$ is the quantum mechanical exchange operator.

Although the stationary condition initially couples the occupied orbitals through the non-diagonal Lagrange multiplier matrix $\boldsymbol{\epsilon}$ with elements $\epsilon_{ji}=\langle\phi_j|\hat{F}|\phi_i\rangle$, we can decouple these equations by transforming the occupied orbitals with a unitary matrix $\mathbf{U}$. Collecting the occupied orbitals into a row vector allows the stationary condition $\hat{F}|\phi_i\rangle=\sum_j^{\text{occ}}\epsilon_{ji}|\phi_j\rangle$ to be written in matrix form as

$$
\hat{F}(|\phi_1\rangle\cdots|\phi_N\rangle)=(|\phi_1\rangle\cdots|\phi_N\rangle)\boldsymbol{\epsilon}
$$

We introduce a new set of transformed occupied orbitals $(|\phi_1'\rangle\cdots|\phi_N'\rangle)$ defined by a unitary transformation matrix $\mathbf{U}$ satisfying $\mathbf{U}^\dagger\mathbf{U}=\mathbf{U}\mathbf{U}^\dagger=\mathbf{I}$ as

$$
(|\phi_1'\rangle\cdots|\phi_N'\rangle)=(|\phi_1\rangle\cdots|\phi_N\rangle)\mathbf{U}
$$

Because $\mathbf{U}$ is unitary, the transformed orbitals remain orthonormal, the physical Slater determinant changes only by an unobservable overall phase factor $\det(\mathbf{U})=e^{i\theta}$, and the density matrix remains completely invariant, $\mathbf{P}'=(\mathbf{C}_{\text{occ}}\mathbf{U})(\mathbf{C}_{\text{occ}}\mathbf{U})^\dagger=\mathbf{C}_{\text{occ}}\mathbf{U}\mathbf{U}^\dagger\mathbf{C}_{\text{occ}}^\dagger=\mathbf{P}$. Since the density matrix is unchanged, the Fock operator $\hat{F}[\mathbf{P}]$ itself is completely unaffected by this rotation. Expressing the original orbitals as $(|\phi_1\rangle\cdots|\phi_N\rangle)=(|\phi_1'\rangle\cdots|\phi_N'\rangle)\mathbf{U}^\dagger$ and substituting them into the coupled operator relation gives

$$
\hat{F}(|\phi_1'\rangle\cdots|\phi_N'\rangle)=(|\phi_1'\rangle\cdots|\phi_N'\rangle)\mathbf{U}^\dagger\boldsymbol{\epsilon}\mathbf{U}
$$

Because the Fock operator is Hermitian, the multiplier matrix is Hermitian ($\epsilon_{ji}=\langle\phi_j|\hat{F}|\phi_i\rangle=\langle\phi_i|\hat{F}|\phi_j\rangle^*=\epsilon_{ij}^*$), and the spectral theorem guarantees that a unitary matrix $\mathbf{U}$ can be chosen whose columns are the eigenvectors of $\boldsymbol{\epsilon}$, diagonalizing it into real orbital energy eigenvalues as $\mathbf{U}^\dagger\boldsymbol{\epsilon}\mathbf{U}=\text{diag}(\epsilon_1,\epsilon_2,\dots,\epsilon_N)$. In this canonical molecular orbital basis, the coupled condition decouples completely into an independent effective one-electron eigenvalue equation for each orbital as

$$
\hat{F}|\phi_i'\rangle=\epsilon_i|\phi_i'\rangle
$$

Dropping the prime notation for simplicity, we substitute the atomic orbital expansion $|\phi_i\rangle=\sum_\nu C_{\nu i}|\chi_\nu\rangle$ into the eigenvalue equation, which gives

$$
\sum_\nu C_{\nu i}\hat{F}|\chi_\nu\rangle=\epsilon_i\sum_\nu C_{\nu i}|\chi_\nu\rangle
$$

Projecting this equation from the left with an atomic orbital basis function $\langle\chi_\mu|$ yields

$$
\sum_\nu C_{\nu i}\langle\chi_\mu|\hat{F}|\chi_\nu\rangle=\epsilon_i\sum_\nu C_{\nu i}\langle\chi_\mu|\chi_\nu\rangle
$$

We identify the matrix elements of the Fock operator in the atomic basis as $F_{\mu\nu}=\langle\chi_\mu|\hat{F}|\chi_\nu\rangle$ and the atomic overlap matrix elements as $S_{\mu\nu}=\langle\chi_\mu|\chi_\nu\rangle$, transforming the projected relation into

$$
\sum_\nu F_{\mu\nu}C_{\nu i}=\epsilon_i\sum_\nu S_{\mu\nu}C_{\nu i}
$$

Evaluating the matrix elements of the Fock operator separates them into the core Hamiltonian and the mean-field potential as

$$
F_{\mu\nu}=H_{\mu\nu}+V_{\mu\nu}^{\text{HF}}=H_{\mu\nu}+J_{\mu\nu}-K_{\mu\nu}
$$

where $H_{\mu\nu}=\langle\chi_\mu|\hat{h}|\chi_\nu\rangle$ is the core Hamiltonian matrix, the local Coulomb matrix elements represent the classical electrostatic repulsion from the electron density,

$$
J_{\mu\nu}=\sum_{\lambda,\sigma}P_{\lambda\sigma}\langle\mu\lambda|\nu\sigma\rangle
$$

and the non-local exchange matrix elements represent the quantum exchange interaction between electrons of identical spin,

$$
K_{\mu\nu}=\sum_{\lambda,\sigma}P_{\lambda\sigma}\langle\mu\lambda|\sigma\nu\rangle
$$

For a single occupied orbital column vector $\mathbf{c}_i$, the algebraic relation is $\mathbf{F}\mathbf{c}_i=\epsilon_i\mathbf{S}\mathbf{c}_i$. Because the Fock matrix $\mathbf{F}$ is completely defined by the occupied density matrix $\mathbf{P}$, diagonalizing $\mathbf{F}$ across the full basis space yields both the $N_{\text{occ}}$ occupied orbitals that minimize the energy and the remaining $N_{\text{AO}}-N_{\text{occ}}$ virtual orbitals spanning the orthogonal complement. Collecting the coefficient vectors for all $N_{\text{AO}}$ orbitals into the columns of matrix $\mathbf{C}$ and the orbital energies $\epsilon_i$ into the diagonal matrix $\mathbf{E}=\text{diag}(\epsilon_1,\epsilon_2,\dots,\epsilon_{N_{\text{AO}}})$ produces the Roothaan–Hall pseudo-eigenvalue equations

$$
\mathbf{F}\mathbf{C}=\mathbf{S}\mathbf{C}\mathbf{E}
$$

which dictate the electronic ground state of the system. Here $\mathbf{F}$ is the Fock matrix, $\mathbf{C}$ is the molecular orbital coefficient matrix, $\mathbf{S}$ is the atomic overlap matrix, and $\mathbf{E}$ is the diagonal matrix of orbital energies. Because the atomic basis functions are not mutually orthogonal, the overlap matrix $\mathbf{S}$ appears on the right-hand side, making this a generalized eigenvalue problem. Furthermore, because the Fock matrix $\mathbf{F}(\mathbf{P})$ depends nonlinearly on the coefficients $\mathbf{C}$ through the density matrix $\mathbf{P}$, the equations cannot be solved in a single step and must instead be solved iteratively using the self-consistent field method until the orbital coefficients and electron density reach numerical convergence. This derivation presents Generalized Hartree–Fock (GHF) in a spin-orbital basis where each state has unit occupancy. For closed-shell systems with spatial orbitals in Restricted Hartree–Fock (RHF), pairs of electrons with opposite spins occupy the same spatial orbital, giving an occupancy factor of two in the density matrix ($P_{\mu\nu}=2\sum_i^{\text{occ}}C_{\mu i}C_{\nu i}^*$) and scaling the exchange term in the Fock matrix by one half because electrons only undergo exchange with electrons of identical spin ($\mathbf{P}_\alpha=\frac{1}{2}\mathbf{P}$), while all electrons contribute to Coulomb repulsion. In the absence of external magnetic fields, the Hamiltonian is time-reversal invariant, so the molecular orbitals can be chosen strictly real, where $\mathbf{U}$ reduces to an orthogonal matrix and $\mathbf{C}^\dagger$ reduces to $\mathbf{C}^T$.

### 7. Canonical Orbitals and Brillouin's Theorem

Once the self-consistent field equations reach convergence, the resulting canonical molecular spin-orbitals diagonalize the Fock operator across the entire orbital space, $\langle\phi_p|\hat{F}|\phi_q\rangle=\epsilon_p\delta_{pq}$. In second quantization, this effective one-electron Fock operator is expressed as

$$
\hat{F}=\sum_p\epsilon_p a_p^\dagger a_p
$$

where $\epsilon_p$ are the canonical orbital energies. Because the Fock operator is strictly diagonal in this canonical basis, any matrix element connecting an occupied spin-orbital $|\phi_i\rangle$ to an unoccupied virtual spin-orbital $|\phi_a\rangle$ vanishes identically, $F_{ia}=\langle\phi_i|\hat{F}|\phi_a\rangle=\epsilon_i\delta_{ia}=0$.

This property leads directly to Brillouin's theorem. In second quantization, any singly excited determinant $|\Phi_i^a\rangle$ is formed by replacing an occupied spin-orbital $i$ with a virtual spin-orbital $a$ as $|\Phi_i^a\rangle=a_a^\dagger a_i|\Phi_0\rangle$. Evaluating the matrix element of the exact Hamiltonian between the Hartree–Fock reference state and this singly excited state using Wick's theorem yields

$$
\langle\Phi_0|\hat{H}|\Phi_i^a\rangle=\langle\Phi_0|\hat{H}a_a^\dagger a_i|\Phi_0\rangle=h_{ia}+\sum_j^{\text{occ}}\left(g_{ijaj}-g_{ijja}\right)=\langle\phi_i|\hat{F}|\phi_a\rangle=F_{ia}=0
$$

which proves that the Hartree–Fock ground state does not interact directly with any singly excited configuration. Brillouin's theorem reflects the variational stationarity of the Hartree–Fock wavefunction, ensuring that small orbital rotations cannot lower the ground-state electronic energy and providing the foundation for post-Hartree–Fock correlation methods.

### 8. Density Fitting for Coulomb and Exchange Potentials (RI-JK)

Evaluating the four-center two-electron Coulomb repulsion matrix elements $J_{\mu\nu}=\sum_{\lambda,\sigma}P_{\lambda\sigma}\langle\mu\lambda|\nu\sigma\rangle$ and exchange elements $K_{\mu\nu}=\sum_{\lambda,\sigma}P_{\lambda\sigma}\langle\mu\lambda|\sigma\nu\rangle$ scales formally as $\mathcal{O}(N_{\text{AO}}^4)$, creating a computational bottleneck for large molecular systems. The resolution-of-the-identity or density fitting approximation expands the product of two primary atomic orbital basis functions in an auxiliary basis set $\{P\}$ as

$$
|\mu\nu\rangle\approx\sum_P C_{\mu\nu}^P|P\rangle
$$

Minimizing the Coulomb metric self-interaction error of the residual density determines the optimal expansion coefficients by solving $\sum_Q J_{PQ}C_{\mu\nu}^Q=\langle\mu P|\nu\rangle$, yielding the fitted product representation

$$
|\mu\nu\rangle\approx\sum_{P,Q}\langle\mu P|\nu\rangle\left[\mathbf{J}^{-1}\right]_{PQ}|Q\rangle
$$

where $J_{PQ}=\iint\frac{P(\mathbf{r}_1)Q(\mathbf{r}_2)}{|\mathbf{r}_1-\mathbf{r}_2|}d\mathbf{r}_1 d\mathbf{r}_2$ is the two-center auxiliary Coulomb metric matrix and $\langle\mu P|\nu\rangle=\iint\frac{\phi_\mu(\mathbf{r}_1)P(\mathbf{r}_2)\phi_\nu(\mathbf{r}_1)}{|\mathbf{r}_1-\mathbf{r}_2|}d\mathbf{r}_1 d\mathbf{r}_2$ denotes the three-center electron repulsion integrals in physicists' Dirac notation. Factorizing the symmetric auxiliary metric as $\mathbf{J}^{-1}=\mathbf{J}^{-1/2}\mathbf{J}^{-1/2}$ allows us to define the transformed three-center intermediate tensor

$$
B_{\mu\nu,Q}=\sum_P\langle\mu P|\nu\rangle\left[\mathbf{J}^{-1/2}\right]_{PQ}
$$

which approximates the four-center two-electron repulsion integrals as $\langle\mu\lambda|\nu\sigma\rangle\approx\sum_Q B_{\mu\nu,Q}B_{\lambda\sigma,Q}$ and reduces the evaluation of the Coulomb matrix to two consecutive matrix-vector products as

$$
q_Q=\sum_{\mu,\nu}B_{\mu\nu,Q}P_{\mu\nu},\quad J_{\mu\nu}=\sum_Q B_{\mu\nu,Q}q_Q
$$

scaling as $\mathcal{O}(N_{\text{AO}}^2 N_{\text{aux}})$ rather than $\mathcal{O}(N_{\text{AO}}^4)$. In Generalized Hartree–Fock, the spatial Coulomb operator is spin-diagonal and contracts with the total density $P_{\mu\nu}^{\text{tot}}=P_{\mu\nu}+P_{\mu+N_{\text{AO}},\nu+N_{\text{AO}}}$, adding the resulting Coulomb potential matrix $\mathbf{J}$ symmetrically to both the $\alpha\alpha$ and $\beta\beta$ spin diagonal blocks of the Fock matrix.

For exact exchange evaluation, the fitting tensor contracts with occupied molecular orbital coefficients to construct the intermediate half-transformed tensor as

$$
W_{\mu i,Q}=\sum_{\nu=1}^{N_{\text{AO}}}C_{\nu i}B_{\mu\nu,Q}
$$

which yields the closed-shell exchange potential matrix through a single symmetric matrix multiplication as

$$
V_{\mu\nu}^{\text{x}}=\frac{1}{2}K_{\mu\nu}=\sum_{i=1}^{N_{\text{occ}}}\sum_Q W_{\mu i,Q}W_{\nu i,Q}=\left[\mathbf{W}\mathbf{W}^T\right]_{\mu\nu}
$$

scaling as $\mathcal{O}(N_{\text{occ}} N_{\text{AO}}^2 N_{\text{aux}})$ and subtracting the exchange potential from the Fock matrix without four-center integral evaluations. In Generalized Hartree–Fock, molecular spin-orbital coefficients partition into spatial spin blocks $\mathbf{C}_{\text{occ}}^\alpha$ and $\mathbf{C}_{\text{occ}}^\beta$, generating intermediate blocks $\mathbf{W}^\alpha$ and $\mathbf{W}^\beta$ whose stacked outer product $\mathbf{K}=\mathbf{W}\mathbf{W}^T$ evaluates all spin blocks of the generalized exchange matrix.

### 9. Density-Fitted Coupled-Perturbed Hartree–Fock Response (RI-CPHF)

Determining first-order orbital responses $\partial\mathbf{C}/\partial X$ and orbital energy derivatives $\partial\boldsymbol{\epsilon}/\partial X$ with respect to a nuclear coordinate $X$ requires solving the Coupled-Perturbed Hartree–Fock (CPHF) linear equations. The orbital rotation matrix $\mathbf{U}^X$ parameterizes the orbital derivative as $\frac{\partial\mathbf{C}}{\partial X}=\mathbf{C}\mathbf{U}^X$, where the occupied–occupied and virtual–virtual blocks are determined directly from the atomic orbital overlap derivative matrix $\mathbf{S}^X=\mathbf{C}^T\frac{\partial\mathbf{S}}{\partial X}\mathbf{C}$ as

$$
U_{ij}^X=-\frac{1}{2}S_{ij}^X,\quad U_{ab}^X=-\frac{1}{2}S_{ab}^X
$$

The occupied–virtual response block $U_{ai}^X$ satisfies the coupled linear system

$$
(\epsilon_a-\epsilon_i)U_{ai}^X+\sum_{b=1}^{N_{\text{vir}}}\sum_{j=1}^{N_{\text{occ}}}A_{ai,bj}U_{bj}^X=-F_{ai}^X+\epsilon_i S_{ai}^X
$$

where $F_{ai}^X=\left[\mathbf{C}^T\frac{\partial\mathbf{F}}{\partial X}\mathbf{C}\right]_{ai}$ is the molecular orbital projection of the unperturbed skeleton Fock derivative matrix and $\mathbf{A}$ represents the electronic Hessian orbital rotation coupling tensor.

In density fitting, all four-center two-electron integral evaluations in both the skeleton Fock derivative and the iterative Hessian-vector contractions are replaced by three-center contractions. The two-electron contribution to the skeleton Fock derivative matrix $\mathbf{F}^{X,\text{2e}}$ is evaluated directly from the derivative fitting tensor $\frac{\partial\mathbf{B}}{\partial X}$ as

$$
F_{\mu\nu}^{X,\text{2e}}=\sum_Q\left(\frac{\partial B_{\mu\nu,Q}}{\partial X}q_Q+B_{\mu\nu,Q}q_Q^X\right)-\sum_{i=1}^{N_{\text{occ}}}\sum_Q\left(\frac{\partial W_{\mu i,Q}}{\partial X}W_{\nu i,Q}+W_{\mu i,Q}\frac{\partial W_{\nu i,Q}}{\partial X}\right)
$$

where $q_Q=\sum_{\mu,\nu}B_{\mu\nu,Q}P_{\mu\nu}$, $q_Q^X=\sum_{\mu,\nu}\frac{\partial B_{\mu\nu,Q}}{\partial X}P_{\mu\nu}$, and $\frac{\partial W_{\mu i,Q}}{\partial X}=\sum_\nu C_{\nu i}\frac{\partial B_{\mu\nu,Q}}{\partial X}$. During each CPHF iteration, the trial perturbation to the occupied orbital coefficients $\mathbf{dC}_x=\mathbf{C}_{\text{vir}}\mathbf{U}_{\text{vo}}^X$ induces the response density $\mathbf{P}^V=2\mathbf{C}_{\text{occ}}\mathbf{dC}_x^T+2\mathbf{dC}_x\mathbf{C}_{\text{occ}}^T$. The response Coulomb potential matrix is evaluated from $\mathbf{P}^V$ via density-fitted Coulomb contraction, while the response exchange potential matrix is evaluated using the perturbed half-transformed tensor $W_{\mu i,Q}^V=\sum_\nu(dC_x)_{\nu i}B_{\mu\nu,Q}$ as

$$
V_{\mu\nu}^{\text{x},V}=\sum_{i=1}^{N_{\text{occ}}}\sum_Q\left(W_{\mu i,Q}^V W_{\nu i,Q}+W_{\mu i,Q}W_{\nu i,Q}^V\right)
$$

scaling as $\mathcal{O}(N_{\text{occ}} N_{\text{AO}}^2 N_{\text{aux}})$ per iteration and eliminating all four-center integrals from orbital response theory.
