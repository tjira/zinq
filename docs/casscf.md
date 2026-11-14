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

---

## IV. Resolution of the Identity Density Fitting

### 7. Density Fitting Representation of Molecular Integrals

To avoid the prohibitive storage and fifth-order scaling of four-index two-electron repulsion integrals, the Resolution of the Identity (RI) or density fitting approximation expands orbital product densities in an auxiliary basis set. The two-electron repulsion integrals in physicist notation are approximated as

$$
\langle pq|rs\rangle=(pr|qs)\approx\sum_P B_{pr}^P B_{qs}^P
$$

where $B_{pq}^P=\sum_Q(pq|Q)[\mathbf{J}^{-1/2}]_{QP}$ is the three-center fitting tensor in the spatial molecular orbital basis, $(pq|Q)$ denotes the three-center Coulomb integral, and $\mathbf{J}$ represents the auxiliary Coulomb metric matrix. In hybrid RI-CASSCF, the molecular orbitals from conventional Hartree–Fock serve as the starting guess, while subsequent macro-iterations utilize the fitting tensor $\mathbf{B}$ directly.

### 8. Generalized Fock Matrix and Configuration Interaction with Density Fitting

Within the density fitting framework, the core Fock matrix is constructed as

$$
F_{pq}^{\text{core}}=h_{pq}+\sum_P B_{pq}^P v_P^{\text{core}}-\sum_{i \in \text{core}}\sum_P B_{pi}^P B_{qi}^P
$$

where $v_P^{\text{core}}=2\sum_{i \in \text{core}}B_{ii}^P$ contracts closed core shells with the auxiliary basis. The core columns of the generalized Fock matrix are evaluated by contracting the active one-particle reduced density matrix as

$$
F_{pj}=2\left(F_{pj}^{\text{core}}+\sum_P B_{pj}^P v_P^{\text{act}}-\frac{1}{2}\sum_{v \in \text{act}}\sum_P B_{p,n_{\text{core}}+v}^P W_{vj}^P\right)
$$

where $v_P^{\text{act}}=\sum_{u,v \in \text{act}}D_{uv}B_{n_{\text{core}}+u,n_{\text{core}}+v}^P$ and $W_{vj}^P=\sum_{u \in \text{act}}D_{uv}B_{n_{\text{core}}+u,j}^P$. The active columns are assembled from the two-particle reduced density matrix as

$$
F_{pu}=\sum_{v \in \text{act}}F_{pv}^{\text{core}}D_{vu}+\sum_{v \in \text{act}}\sum_P B_{p,n_{\text{core}}+v}^P Q_{uv}^P
$$

where $Q_{uv}^P=\sum_{w,x \in \text{act}}d_{uvwx}B_{n_{\text{core}}+w,n_{\text{core}}+x}^P$. Active-space configuration interaction integrals are evaluated on the fly by contracting the transformed three-center tensors, completely bypassing the formation of four-index electron repulsion tensors.

---

## V. Analytic Nuclear Gradients with Density Fitting

### 9. Dual-Number Forward Automatic Differentiation

The analytic nuclear gradient of the converged single-state CASSCF energy with respect to a nuclear coordinate $x$ is evaluated by differentiating the electronic energy functional while enforcing the molecular orbital orthonormality constraints. Because the CASSCF energy is fully variational with respect to both the configuration interaction coefficients $\mathbf{c}$ and non-redundant orbital rotations $\boldsymbol{\kappa}$ at convergence, the Hellmann–Feynman theorem guarantees that the implicit response contributions $\frac{\partial E}{\partial\mathbf{c}}\frac{\partial\mathbf{c}}{\partial x}=0$ and $\frac{\partial E}{\partial\boldsymbol{\kappa}}\frac{\partial\boldsymbol{\kappa}}{\partial x}=0$ vanish identically, eliminating the need to solve coupled-perturbed multiconfigurational equations.

The geometric derivatives of the molecular orbital coefficients originate solely from preserving the orbital orthonormality condition $\mathbf{C}^T\mathbf{S}\mathbf{C}=\mathbf{I}$, yielding the connection matrix derivative

$$
\frac{\partial\mathbf{C}}{\partial x}=-\frac{1}{2}\mathbf{C}\left(\mathbf{C}^T\frac{\partial\mathbf{S}}{\partial x}\mathbf{C}\right)
$$

which is assembled using the derivative atomic orbital overlap matrix $\frac{\partial\mathbf{S}}{\partial x}$. Forward-mode automatic differentiation represents each atomic orbital quantity as a dual number $\tilde{a}=a+\epsilon\frac{\partial a}{\partial x}$ with dual unit $\epsilon^2=0$, bundling the baseline integrals and nuclear derivatives into unified algebraic objects.

### 10. Density Fitting Three-Center Tensor Derivatives

In Resolution of the Identity density fitting, the three-center fitting tensor in the atomic orbital basis is given by $\mathbf{B}=\mathbf{I}_{3c}\mathbf{J}^{-1/2}$, where $\mathbf{I}_{3c}$ contains the three-center two-electron repulsion integrals and $\mathbf{J}$ represents the Coulomb metric matrix. Differentiating this fitting tensor with respect to nuclear coordinate $x$ via the product rule yields

$$
\frac{\partial B_{\mu\nu}^P}{\partial x}=\sum_Q\frac{\partial(\mu\nu|Q)}{\partial x}[\mathbf{J}^{-1/2}]_{QP}+\sum_Q(\mu\nu|Q)\frac{\partial[\mathbf{J}^{-1/2}]_{QP}}{\partial x}
$$

where the derivative of the inverse square-root Coulomb metric $\frac{\partial\mathbf{J}^{-1/2}}{\partial x}$ is evaluated from the spectral decomposition of $\mathbf{J}=\mathbf{U}\mathbf{W}\mathbf{U}^T$ and derivative metric integrals $\frac{\partial\mathbf{J}}{\partial x}$.

The atomic orbital dual core Hamiltonian $\tilde{\mathbf{H}}$ and dual fitting tensor $\tilde{\mathbf{B}}$ are subsequently transformed to the molecular orbital basis using the dual coefficient matrix $\tilde{\mathbf{C}}=\mathbf{C}+\epsilon\frac{\partial\mathbf{C}}{\partial x}$ through forward dual integral transformations. The active-space configuration interaction Hamiltonian elements between Slater determinants $|\Phi_m\rangle$ and $|\Phi_n\rangle$ are evaluated with dual Slater–Condon rules, directly yielding the total nuclear gradient

$$
\frac{\partial E}{\partial x}=\frac{\partial V_{\text{nuc}}}{\partial x}+\sum_{m \le n}(2-\delta_{mn})c_m c_n\left[\frac{\partial\langle\Phi_m|\hat{H}|\Phi_n\rangle}{\partial x}\right]_{\epsilon}
$$

where $\frac{\partial V_{\text{nuc}}}{\partial x}$ is the nuclear repulsion gradient and $[\cdot]_\epsilon$ extracts the first-order dual derivative component.

### 11. State-Averaged Excited State Gradients and Coupled-Perturbed MCSCF

For state-averaged CASSCF with multiple target electronic states, the molecular orbitals are optimized for the weighted ensemble average energy $\bar{E}=\sum_I w_I E_I$, meaning that the orbital gradient $\mathbf{g}^{(s)}$ of an individual electronic state $s$ does not vanish. To incorporate the resulting orbital relaxation without explicitly evaluating derivative orbital coefficients, the coupled-perturbed MCSCF Z-vector response equation is solved as

$$
\mathbf{A}_{\text{SA}}\mathbf{Z}^{(s)}=-\mathbf{g}^{(s)}
$$

where $\mathbf{A}_{\text{SA}}$ represents the orbital Hessian of the state-averaged energy and $\mathbf{Z}^{(s)}$ is the orbital response Lagrange multiplier vector.

The total analytical gradient for target state $s$ is obtained by appending the orbital relaxation contraction to the dual Hellmann–Feynman derivative as

$$
\frac{\partial E_s}{\partial x}=\left[\frac{\partial E_s}{\partial x}\right]_{\text{HF}}-\sum_k Z_k^{(s)}B_k^{[x]}
$$

where $B_k^{[x]}=2(dF_{pq}-dF_{qp})$ is the derivative generalized Fock response assembled from the dual core Hamiltonian derivative $d\mathbf{h}_{\text{mo}}$ and derivative density fitting two-electron integrals $dg_{pqrs}=\sum_P\left(dB_{pr}^P B_{qs}^P+B_{pr}^P dB_{qs}^P\right)$.
