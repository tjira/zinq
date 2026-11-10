# Complete Active Space Self-Consistent Field

This document provides a clear, simple, and mathematically rigorous explanation of the Complete Active Space Self-Consistent Field (CASSCF) method in second quantization. Designed to treat strong static electron correlation, near-degeneracies, and open-shell systems, this guide explains how the orbital space is partitioned, how the multiconfigurational wavefunction and reduced density matrices are formulated, how the generalized Fock matrix and orbital gradient are derived, and how two-step orbital and configuration optimization is achieved.

---

## I. Orbital Space Partitioning and Multiconfigurational Wavefunction

### 1. Partitioning of the Molecular Orbital Space

Single-determinant methods like Hartree–Fock break down fundamentally for molecules with near-degenerate electronic states, broken bonds, or open-shell diradical character because no single Slater determinant can qualitatively describe the true wavefunction. The Complete Active Space Self-Consistent Field (CASSCF) method resolves this static correlation problem by partitioning the $M$ spatial molecular orbitals into three disjoint physical subspaces according to their occupation characteristics.

The inactive or core subspace contains $n_{\text{core}}$ spatial orbitals that remain doubly occupied across all electronic configurations; the active subspace contains $n_{\text{act}}$ spatial orbitals accommodating $n_{\text{act,el}}$ active electrons whose occupations vary freely between zero, one, and two; and the secondary or virtual subspace contains $n_{\text{virt}}$ unoccupied spatial orbitals that remain empty in all configurations. The total number of electrons in the molecule is partitioned as $N_{\text{elec}}=2n_{\text{core}}+n_{\text{act,el}}$, while the total number of molecular orbitals satisfies $M=n_{\text{core}}+n_{\text{act}}+n_{\text{virt}}$.

### 2. Complete Active Space Configuration Expansion

Within the active orbital subspace, the electronic wavefunction $|\Psi_k\rangle$ for electronic state $k$ is constructed as a Full Configuration Interaction (FCI) expansion over all symmetry- and spin-adapted active Slater determinants as

$$
|\Psi_k\rangle=\sum_a c_{ak}|\Phi_a\rangle
$$

where $c_{ak}$ are the state Configuration Interaction coefficients satisfying the normalization condition $\sum_a c_{ak}^2=1$. In second quantization, each determinant $|\Phi_a\rangle$ is formed by appending the fully occupied core spin-orbitals to an active spin-orbital creation string acting on the vacuum state as

$$
|\Phi_a\rangle=\left(\prod_{i \in \text{core}} a_{i\alpha}^\dagger a_{i\beta}^\dagger\right)\left(\prod_{u \in \text{act}_a} a_{u\sigma}^\dagger\right)|0\rangle
$$

By allowing the active electrons to distribute across all active orbitals in every possible spin-conserving combination, CASSCF captures the complete static electron correlation within the selected valence active space.

---

## II. Reduced Density Matrices and Electronic Energy

### 3. Spin-Summed Excitation Operators and Reduced Density Matrices

To express the electronic energy and orbital gradients independently of specific spin projections, we employ the spin-summed one-body and two-body excitation operators defined as

$$
\hat{E}_{pq}=\sum_\sigma a_{p\sigma}^\dagger a_{q\sigma},\quad\hat{e}_{pqrs}=\sum_{\sigma,\tau}a_{p\sigma}^\dagger a_{r\tau}^\dagger a_{s\tau} a_{q\sigma}=\hat{E}_{pq}\hat{E}_{rs}-\delta_{qr}\hat{E}_{ps}
$$

The spatial distributions of electrons in the active space are characterized by the active one-particle reduced density matrix (1-RDM) $\mathbf{D}$ and two-particle reduced density matrix (2-RDM) $\mathbf{d}$. For state $k$, their elements are defined by the active space expectation values

$$
D_{uv}^{(k)}=\langle\Psi_k|\hat{E}_{uv}|\Psi_k\rangle=\sum_{a,b}c_{ak}c_{bk}\langle\Phi_a|\hat{E}_{uv}|\Phi_b\rangle
$$

and

$$
d_{uvwx}^{(k)}=\langle\Psi_k|\hat{e}_{uvwx}|\Psi_k\rangle=\sum_{a,b}c_{ak}c_{bk}\langle\Phi_a|\hat{e}_{uvwx}|\Phi_b\rangle
$$

In State-Averaged CASSCF (SA-CASSCF), the molecular orbitals are optimized for an ensemble of $M_{\text{states}}$ target electronic states using normalized non-negative state weights $w_I$ satisfying $\sum_I w_I=1$. The state-averaged active reduced density matrices are constructed as the weighted statistical averages

$$
\mathbf{D}=\sum_{I=1}^{M_{\text{states}}}w_I\mathbf{D}^{(I)},\quad\mathbf{d}=\sum_{I=1}^{M_{\text{states}}}w_I\mathbf{d}^{(I)}
$$

which prevents root-flipping instabilities during numerical optimization and yields a single balanced, orthogonal orbital set capable of describing ground and excited states on an equal footing. For single-state CASSCF, $M_{\text{states}}=1$ and $w_1=1$, recovering the state-specific density matrices $\mathbf{D}=\mathbf{D}^{(1)}$ and $\mathbf{d}=\mathbf{d}^{(1)}$ directly.

### 4. Electronic Energy Functional in the Molecular Basis

Transforming the electronic Hamiltonian into the spatial molecular orbital basis using the spin-summed excitation operators gives $\hat{H}=\sum_{pq}h_{pq}\hat{E}_{pq}+\frac{1}{2}\sum_{pqrs}\langle pq|rs\rangle\sum_{\sigma,\tau}a_{p\sigma}^\dagger a_{q\tau}^\dagger a_{s\tau}a_{r\sigma}$, where $\langle pq|rs\rangle=\int\int\phi_p^*(\mathbf{r}_1)\phi_q^*(\mathbf{r}_2)\frac{1}{r_{12}}\phi_r(\mathbf{r}_1)\phi_s(\mathbf{r}_2)d\mathbf{r}_1d\mathbf{r}_2$ denotes the spatial two-electron repulsion integral in physicist notation. Evaluating the expectation value of this Hamiltonian over the multiconfigurational wavefunction isolates the inactive core contributions, expressing the total electronic energy as

$$
E=E_{\text{core}}+\sum_{u,v \in \text{act}}D_{uv}F_{uv}^{\text{core}}+\frac{1}{2}\sum_{u,v,w,x \in \text{act}}d_{uvwx}\langle uw|vx\rangle
$$

In this energy expression, $F_{pq}^{\text{core}}$ is the effective one-electron core Fock Hamiltonian incorporating the mean-field Coulomb and exchange interactions from the doubly occupied inactive electrons, defined in physicist notation as

$$
F_{pq}^{\text{core}}=h_{pq}+\sum_{i \in \text{core}}\left(2\langle pi|qi\rangle-\langle pi|iq\rangle\right)
$$

where $h_{pq}=\langle\phi_p|\hat{h}|\phi_q\rangle$ is the core Hamiltonian from Hartree–Fock, and $E_{\text{core}}=\sum_{i \in \text{core}}\left(h_{ii}+F_{ii}^{\text{core}}\right)$ represents the energy contribution of the closed inactive shell.

---

## III. Generalized Brillouin Theorem and Orbital Optimization

### 5. Generalized Fock Matrix and Orbital Gradient

Unlike Hartree–Fock theory where only molecular orbitals are optimized, or Configuration Interaction where only CI expansion coefficients are optimized, CASSCF requires simultaneous variational optimization of both the configuration coefficients $\mathbf{c}_k$ and the molecular orbital coefficients $\mathbf{C}$. Variations in the molecular orbitals are generated by unitary transformations $\mathbf{C}(\boldsymbol{\kappa})=\mathbf{C}e^{\boldsymbol{\kappa}}$, parameterized by an anti-symmetric matrix of orbital rotation parameters $\kappa_{pq}=-\kappa_{qp}$.

Differentiating the electronic energy with respect to an orbital rotation parameter $\kappa_{pq}$ defines the orbital gradient matrix $\mathbf{g}$ as

$$
g_{pq}=2(F_{pq}-F_{qp})
$$

where $\mathbf{F}$ is the generalized Fock matrix. The elements of the generalized Fock matrix are evaluated block-wise according to the subspace of the column orbital index in physicist notation as

$$
F_{pj}=2F_{pj}^{\text{core}}+2\sum_{u,v \in \text{act}}D_{uv}\left(\langle pu|jv\rangle-\frac{1}{2}\langle pu|vj\rangle\right)
$$

for inactive core columns $j \in \text{core}$, and as

$$
F_{pu}=\sum_{v \in \text{act}}F_{pv}^{\text{core}}D_{vu}+\sum_{v,w,x \in \text{act}}d_{uvwx}\langle pw|vx\rangle
$$

for active columns $u \in \text{act}$, while vanishing identically for virtual columns $a \in \text{virt}$ as

$$
F_{pa}=0
$$

At stationary points of the energy hypersurface, the generalized Brillouin theorem requires that the generalized Fock matrix is symmetric ($F_{pq}=F_{qp}$), which causes the orbital gradient $g_{pq}$ to vanish for all non-redundant orbital pairs. Rotations between orbitals within the same subspace (core–core, active–active, and virtual–virtual) leave the physical CASSCF energy strictly invariant and have zero gradient, meaning that non-redundant rotations occur exclusively between different orbital subspaces (core–active, core–virtual, and active–virtual).

### 6. Macro-Iteration Algorithm and Unitary Cayley Transform

The CASSCF optimization proceeds via a two-step macro-iteration cycle that alternates between configuration interaction and orbital optimization until self-consistency is reached. In the CI step, the molecular orbitals are held fixed, and the active-space configuration interaction Hamiltonian is diagonalized to determine state energies $E_k$, expansion eigenvectors $\mathbf{c}_k$, and active reduced density matrices $\mathbf{D}$ and $\mathbf{d}$.

In the orbital step, the density matrices are held fixed, and the generalized Fock matrix $\mathbf{F}$ and orbital gradient $\mathbf{g}$ are assembled. An approximate quasi-Newton step computes the non-redundant orbital rotation parameters $\boldsymbol{\kappa}$ using a diagonal Hessian approximation with diagonal orbital energies $F_{pp}$ as

$$
\kappa_{bj}=-\frac{g_{bj}}{4(F_{bb}-F_{jj})},\quad\kappa_{bu}=-\frac{g_{bu}}{2D_{uu}(F_{bb}-F_{uu})},\quad\kappa_{uj}=-\frac{g_{uj}}{2(2-D_{uu})(F_{uu}-F_{jj})}
$$

To preserve the strict orthonormality constraint $\mathbf{C}^\dagger\mathbf{S}\mathbf{C}=\mathbf{I}$ without requiring expensive matrix exponentials, the anti-symmetric generator matrix $\boldsymbol{\kappa}$ is mapped to an orthogonal orbital rotation matrix $\mathbf{U}$ via the Cayley transform as

$$
\mathbf{U}=\left(\mathbf{I}-\frac{1}{2}\boldsymbol{\kappa}\right)^{-1}\left(\mathbf{I}+\frac{1}{2}\boldsymbol{\kappa}\right)
$$

which is solved efficiently by LU factorization. In the absence of external magnetic fields, the molecular orbitals can be chosen strictly real, where $\mathbf{C}^\dagger\mathbf{S}\mathbf{C}=\mathbf{I}$ reduces to $\mathbf{C}^T\mathbf{S}\mathbf{C}=\mathbf{I}$ and $\mathbf{U}$ is an orthogonal matrix. The molecular orbital coefficients are updated as $\mathbf{C}\leftarrow\mathbf{C}\mathbf{U}$, and the macro-iteration cycle is repeated until the maximum orbital gradient element $\max|g_{pq}|$ and the energy change drop below the convergence threshold, completing the self-consistent multiconfigurational solution.
