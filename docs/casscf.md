# Complete Active Space Self-Consistent Field

This document provides a detailed and mathematically rigorous explanation of the Complete Active Space Self-Consistent Field (CASSCF) method and its computational implementation in Zinq. The guide presents the multiconfigurational wavefunction formulation, two-step orbital and configuration interaction optimization, generalized Fock matrix construction, Cayley unitary orbital rotation, and analytic nuclear derivatives.

---

## I. Orbital Space Partitioning and Wavefunction Formulation

### 1. Partitioning of the Molecular Orbital Space

In molecular systems exhibiting near-degeneracy and strong static electron correlation, such as bond dissociation processes, diradicals, and open-shell transition metal complexes, a single Slater determinant cannot qualitatively describe the electronic structure. The Complete Active Space Self-Consistent Field method addresses this multiconfigurational problem by partitioning the molecular orbital space into three distinct subspaces according to physical occupation. The inactive or core space contains $n_{\text{core}}$ spatial orbitals that remain doubly occupied across all electronic configurations; the active space contains $n_{\text{act}}$ spatial orbitals accommodating $n_{\text{act,el}}$ active electrons whose occupations vary freely between zero, one, and two; and the virtual or secondary space contains $n_{\text{virt}}$ unoccupied spatial orbitals that remain empty in every determinant.

### 2. Complete Active Space Configuration Expansion

Within the active orbital subspace, the electronic wavefunction $|\Psi_k\rangle$ for state $k$ is constructed as an exact Full Configuration Interaction (FCI) expansion over all symmetry- and spin-adapted Slater determinants $|\Phi_a\rangle$ formed by distributing the $n_{\text{act,el}}$ active electrons among the $n_{\text{act}}$ active orbitals as

$$
|\Psi_k\rangle=\sum_a C_{ak}|\Phi_a\rangle
$$

where $C_{ak}$ are the state configuration interaction expansion coefficients satisfying the normalization condition $\sum_a C_{ak}^2=1$. The total number of electrons in the molecular system is partitioned as $N_{\text{elec}}=2n_{\text{core}}+n_{\text{act,el}}$, and each determinant $|\Phi_a\rangle$ is formed by appending the fully occupied core spin-orbitals to the active spin-orbital excitation string, ensuring that core electrons are preserved while all active configurations match the spin multiplicity $2S+1$.

---

## II. Electronic Hamiltonian and Reduced Density Matrices

### 3. Active Space Integrals and CI Diagonalization

To determine the state CI expansion coefficients $\mathbf{C}_{\text{ci}}$ and electronic energies for the current molecular orbital basis $\mathbf{C}$, the one-electron core Hamiltonian $\mathbf{h}$ and two-electron Coulomb repulsion integrals $(\mu\nu|\lambda\sigma)$ are transformed into the molecular orbital basis. The effective one-electron core Fock Hamiltonian incorporating the mean-field Coulomb and exchange interactions of the doubly occupied inactive core electrons is evaluated as

$$
F_{pq}^{\text{core}}=h_{pq}+\sum_{i \in \text{core}}\left[2(pq|ii)-(pi|iq)\right]
$$

where $h_{pq}=\sum_{\mu\nu}C_{\mu p}C_{\nu q}h_{\mu\nu}$, and $(pq|rs)=\sum_{\mu\nu\lambda\sigma}C_{\mu p}C_{\nu q}C_{\lambda r}C_{\sigma s}(\mu\nu|\lambda\sigma)$ in chemist notation. In the basis of active determinants, the matrix elements of the electronic Hamiltonian operator $\hat{H}$ are evaluated using the Slater–Condon rules with the core Fock operator and active molecular integrals. Diagonalizing the active configuration interaction Hamiltonian matrix via symmetric eigensolvers yields the state energies $E_k$ and the corresponding eigenvector columns $\mathbf{C}_{\text{ci}}$.

### 4. Active One- and Two-Particle Reduced Density Matrices

Orbital optimization requires the spin-summed one-particle reduced density matrix $\mathbf{D}$ and two-particle reduced density matrix $\mathbf{d}$ defined in the active orbital space. Using the spin-summed excitation operators $\hat{E}_{uv}=\sum_\sigma a_{u\sigma}^\dagger a_{v\sigma}$ and $\hat{e}_{uvwx}=\sum_{\sigma\tau}a_{u\sigma}^\dagger a_{w\tau}^\dagger a_{x\tau}a_{v\sigma}=\hat{E}_{uv}\hat{E}_{wx}-\delta_{vw}\hat{E}_{ux}$, the elements of the active one-particle reduced density matrix for electronic state $k$ are defined as

$$
D_{uv}=\langle\Psi_k|\hat{E}_{uv}|\Psi_k\rangle=\sum_{a,b}C_{ak}C_{bk}\langle\Phi_a|\hat{E}_{uv}|\Phi_b\rangle
$$

where the transition elements are non-zero only when determinants $a$ and $b$ differ by at most one active spin-orbital excitation. The elements of the active two-particle reduced density matrix are evaluated as

$$
d_{uvwx}=\langle\Psi_k|\hat{e}_{uvwx}|\Psi_k\rangle=\sum_{a,b}C_{ak}C_{bk}\langle\Phi_a|\hat{e}_{uvwx}|\Phi_b\rangle
$$

which quantify correlated two-electron pair distributions among active orbitals. The total electronic energy of the multiconfigurational state is expressed in terms of these reduced density matrices as

$$
E_{\text{elec}}=E_{\text{core}}+\sum_{uv \in \text{act}}D_{uv}F_{uv}^{\text{core}}+\frac{1}{2}\sum_{uvwx \in \text{act}}d_{uvwx}(uv|wx)
$$

where $E_{\text{core}}=\sum_{i \in \text{core}}\left(h_{ii}+F_{ii}^{\text{core}}\right)$ represents the energy contribution of the doubly occupied inactive core electrons.

---

## III. Generalized Fock Matrix and Orbital Gradient

### 5. Generalized Fock Matrix Construction

In standard Hartree–Fock theory, the Brillouin theorem states that singly excited determinants do not interact with the Hartree–Fock reference determinant, which is equivalent to requiring the Fock operator to be diagonal. In multiconfigurational self-consistent field theory, this condition is generalized through the generalized Fock matrix $\mathbf{F}^{\text{gen}}$. In conventional MCSCF theory, the elements of the generalized Fock matrix $F_{pq}^{\text{gen}}$ are classified according to the subspace of the column orbital index $q$ as

$$
F_{pj}^{\text{gen}}=2F_{pj}^{\text{tot}}=2\left(F_{pj}^{\text{core}}+\sum_{uv \in \text{act}}D_{uv}\left[(pj|uv)-\frac{1}{2}(pv|uj)\right]\right)
$$

for inactive core columns $j \in \text{core}$, where $F_{pj}^{\text{tot}}$ combines the inactive core potential with the Coulomb and exchange fields generated by the active electron density. For active columns $u \in \text{act}$, the generalized Fock matrix elements are computed from the active one-particle reduced density matrix contracted with the core Fock operator and the active two-particle reduced density matrix contracted with the two-electron integrals as

$$
F_{pu}^{\text{gen}}=\sum_{v \in \text{act}}F_{pv}^{\text{core}}D_{vu}+\sum_{vwx \in \text{act}}d_{uvwx}(pv|wx)
$$

which represents the effective mean field experienced by an active electron interacting with the correlated multiconfigurational active space. For secondary virtual columns $a \in \text{virt}$, the generalized Fock matrix elements vanish identically as

$$
F_{pa}^{\text{gen}}=0
$$

because secondary virtual orbitals are unpopulated in every Slater determinant of the wavefunction.

### 6. Orbital Gradient and the Generalized Brillouin Theorem

The variation of the multiconfigurational electronic energy with respect to an infinitesimal anti-symmetric orbital rotation parameter $\kappa_{pq}=-\kappa_{qp}$ defines the orbital gradient matrix $\mathbf{g}$ as

$$
g_{pq}=2\left(F_{pq}^{\text{gen}}-F_{qp}^{\text{gen}}\right)
$$

which provides the driving force for molecular orbital optimization. At stationary points of the CASSCF energy hypersurface, the generalized Brillouin theorem requires that the generalized Fock matrix is symmetric ($F_{pq}^{\text{gen}}=F_{qp}^{\text{gen}}$), causing the orbital gradient $g_{pq}$ to vanish for all non-redundant orbital pairs. Non-redundant orbital rotations occur exclusively between different orbital subspaces, specifically core–active, core–virtual, and active–virtual orbital pairs, whereas rotations within the same subspace (core–core, active–active, and virtual–virtual) leave the physical CASSCF energy invariant and have zero orbital gradients.

---

## IV. Orbital Optimization and Macro-Iteration Algorithm

### 7. Approximate Quasi-Newton Step and Unitary Cayley Transform

To update the molecular orbitals toward self-consistency without constructing the exact four-index orbital Hessian matrix, an approximate diagonal quasi-Newton step is evaluated. The non-redundant orbital rotation parameters $\kappa_{pq}$ are computed using the diagonal Hessian approximation as

$$
\kappa_{bj}=-\frac{g_{bj}}{4(F_{bb}^{\text{diag}}-F_{jj}^{\text{diag}})}
$$

for core–virtual rotations between core orbital $j$ and virtual orbital $b$, and as

$$
\kappa_{bu}=-\frac{g_{bu}}{2D_{uu}(F_{bb}^{\text{diag}}-F_{uu}^{\text{diag}})}
$$

for active–virtual rotations between active orbital $u$ and virtual orbital $b$, and as

$$
\kappa_{uj}=-\frac{g_{uj}}{2(2-D_{uu})(F_{uu}^{\text{diag}}-F_{jj}^{\text{diag}})}
$$

for core–active rotations between core orbital $j$ and active orbital $u$, where $D_{uu}$ is the diagonal occupation of active orbital $u$. To preserve the strict orthonormality condition $\mathbf{C}^T\mathbf{S}\mathbf{C}=\mathbf{I}$ throughout optimization, the anti-symmetric generator matrix $\boldsymbol{\kappa}$ is mapped to an orthogonal rotation matrix $\mathbf{U}$ via the Cayley transform defined as

$$
\mathbf{U}=\left(\mathbf{I}-\frac{1}{2}\boldsymbol{\kappa}\right)^{-1}\left(\mathbf{I}+\frac{1}{2}\boldsymbol{\kappa}\right)
$$

which is solved numerically by LU factorization of the linear system $\left(\mathbf{I}-\frac{1}{2}\boldsymbol{\kappa}\right)\mathbf{U}=\mathbf{I}+\frac{1}{2}\boldsymbol{\kappa}$. The molecular orbital coefficient matrix is updated as

$$
\mathbf{C}^{(m+1)}=\mathbf{C}^{(m)}\mathbf{U}
$$

and the macro-iteration cycle of integral transformation, CAS-CI diagonalization, reduced density matrix evaluation, generalized Fock construction, and Cayley orbital rotation is repeated until the maximum orbital gradient element $\max|g_{pq}|$ and the energy change drop below the convergence threshold.

### 8. Orbital Canonicalization

Because the CASSCF energy is invariant under arbitrary unitary rotations within the core, active, and virtual subspaces, the converged molecular orbitals within each subspace are not uniquely determined. To obtain canonical, symmetry-adapted molecular orbitals, the core–core and virtual–virtual diagonal blocks of the core Fock operator are diagonalized independently after macro-iteration convergence. This transformation orders the molecular orbitals within each subspace according to increasing orbital energy eigenvalues while leaving the total multiconfigurational wavefunction and energy completely invariant.

---

## V. Analytic Nuclear Gradients and Nonadiabatic Couplings

### 9. Automatic Differentiation and the Generalized Hellmann–Feynman Theorem

Analytic nuclear gradients of the CASSCF state energy with respect to Cartesian nuclear coordinates $\mathbf{R}_A$ are evaluated using forward-mode automatic differentiation with dual numbers. By the generalized Hellmann–Feynman theorem, because the converged CASSCF wavefunction is stationary with respect to all configuration interaction coefficients and non-redundant orbital rotations, explicit coupled-perturbed multiconfigurational orbital response contributions vanish from the energy gradient. The nuclear gradient is expressed in conventional form as

$$
\frac{\partial E_k}{\partial R_{A\alpha}}=\sum_{pq}D_{pq}\frac{\partial h_{pq}}{\partial R_{A\alpha}}+\frac{1}{2}\sum_{pqrs}d_{pqrs}\frac{\partial(pq|rs)}{\partial R_{A\alpha}}-\sum_{pq}F_{pq}^{\text{gen}}\left(\mathbf{C}^T\frac{\partial\mathbf{S}}{\partial R_{A\alpha}}\mathbf{C}\right)_{pq}+\frac{\partial V_{\text{NN}}}{\partial R_{A\alpha}}
$$

where $V_{\text{NN}}$ is the classical nuclear repulsion energy. In Zinq, this gradient is evaluated exactly to machine precision by lifting atomic orbital integrals to the dual number algebra with the orbital connection derivative $d\mathbf{C}=-\frac{1}{2}\mathbf{C}(\mathbf{C}^T d\mathbf{S}\,\mathbf{C})$ maintaining basis set orthonormality.

### 10. Nonadiabatic Coupling Vectors

For nonadiabatic molecular dynamics and electronic transitions involving multiple states, the first-derivative nonadiabatic coupling vector between states $J$ and $K$ describes the electronic momentum coupling defined as

$$
\mathbf{d}_{JK}(\mathbf{R})=\langle\Psi_J|\boldsymbol{\nabla}_{\mathbf{R}}\Psi_K\rangle=\frac{\langle\Psi_J|(\boldsymbol{\nabla}_{\mathbf{R}}\hat{H})|\Psi_K\rangle}{E_K-E_J}
$$

which governs radiationless transitions between potential energy surfaces. In the CASSCF framework, nonadiabatic coupling vectors are evaluated numerically via central finite difference coordinate displacements as

$$
\mathbf{d}_{JK}(\mathbf{R})=\lim_{h \to 0}\frac{\langle\Psi_J(\mathbf{R}-h\hat{\mathbf{e}})|\Psi_K(\mathbf{R}+h\hat{\mathbf{e}})\rangle-\langle\Psi_J(\mathbf{R}+h\hat{\mathbf{e}})|\Psi_K(\mathbf{R}-h\hat{\mathbf{e}})\rangle}{4h}
$$

combined with active determinant overlap evaluation and molecular orbital phase tracking to preserve wavefunction continuity for Fewest Switches Surface Hopping and Ehrenfest dynamics.