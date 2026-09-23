# Population Analysis

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of population analysis and bond order calculations as implemented in our scientific computing framework. If you want to understand how we can divide a continuous quantum mechanical electron cloud among the individual atoms in a molecule to calculate atomic charges and evaluate covalent bond orders between pairs of atoms, this guide is written step-by-step for you.

---

## I. Mulliken Population Analysis

To understand population analysis, we must first look at how quantum mechanics describes electrons in a molecule. Unlike classical physics, where electrons are treated as tiny charged balls with definite positions, quantum mechanics describes electrons as a continuous cloud of probability density spread out over the entire molecule. However, chemists prefer to think of molecules as collections of individual atoms that share electrons through chemical bonds and carry partial electrical charges. Population analysis is a mathematical method used to divide this continuous electron density cloud and assign a specific number of electrons to each individual atom. Mulliken population analysis is the most common method to achieve this.

### 1. Basis Function and Atomic Populations

In quantum chemistry, we build the molecular orbitals using atomic orbital basis functions. Because these basis functions are centered on different atoms and extend out into space, they overlap with each other. This means that a portion of the electron density lives in the space between the atoms and is shared between them. Mulliken population analysis divides the electronic charge by taking this shared overlap density and splitting it exactly in half between the two sharing atoms. The net electron population $N_\mu$ associated with a single basis function $\chi_\mu$ is calculated using the elements of the density matrix $\mathbf{P}$ and the overlap matrix $\mathbf{S}$ as

$$
N_\mu=\sum_\nu P_{\mu\nu}S_{\mu\nu}
$$

which mathematically corresponds to a diagonal element of the matrix product $\mathbf{P}\mathbf{S}$. If we were to multiply the entire density matrix $\mathbf{P}$ by the overlap matrix $\mathbf{S}$ just to extract the diagonal elements, the calculation would be slow and require a lot of memory. To optimize this, our codebase implements a trick where it directly computes the dot product of the corresponding rows of the density and overlap matrices, which calculates the diagonal elements in linear time. The total electron population $N_A$ assigned to atom $A$ is then obtained by summing the populations of all the basis functions that belong to that atom as

$$
N_A=\sum_{\mu\in A}N_\mu
$$

where the code maps each basis function to its corresponding atomic center using a pre-constructed lookup array named `sys.bf2at` representing the basis function to atom association.

### 2. Net Atomic Charges

Once we have calculated the total electron population $N_A$ for each atom, we can find the net atomic charge $q_A$ of atom $A$. The net charge is calculated by subtracting the total electron population from the nuclear charge $Z_A$, which is the number of protons in the nucleus, using the equation

$$
q_A=Z_A-N_A
$$

where a positive net charge $q_A$ indicates that the atom has lost electron density and carries a partial positive charge, while a negative net charge indicates that the atom has pulled electron density towards itself and carries a partial negative charge. These charges are useful for understanding the polarity of the molecule, its dipole moment, and how it will interact with other molecules.

---

## II. Löwdin Population Analysis

While Mulliken population analysis provides a simple way to partition the electron density, it can suffer from unphysical charge distributions when large or diffuse basis sets are used because it partitions off-diagonal overlap elements equally without orthogonalization. Löwdin population analysis overcomes this limitation by transforming the non-orthogonal atomic orbital basis into a symmetrically orthogonalized basis before assigning electronic populations.

### 1. Symmetric Orthogonalization and Populations

Löwdin symmetric orthogonalization constructs an orthogonal basis using the square root of the overlap matrix $\mathbf{S}^{1/2}$. The symmetric square root matrix $\mathbf{S}^{1/2}$ is computed by diagonalizing the overlap matrix $\mathbf{S}=\mathbf{U}\boldsymbol{\Lambda}\mathbf{U}^T$ through eigenvalue decomposition and evaluating

$$
\mathbf{S}^{1/2}=\mathbf{U}\boldsymbol{\Lambda}^{1/2}\mathbf{U}^T
$$

where $\boldsymbol{\Lambda}^{1/2}=\text{diag}(\sqrt{\lambda_1},\dots,\sqrt{\lambda_N})$ contains the square roots of the overlap eigenvalues. The electron density matrix expressed in this symmetrically orthogonalized Löwdin basis is given by

$$
\mathbf{P}^{\text{Löwdin}}=\mathbf{S}^{1/2}\mathbf{P}\mathbf{S}^{1/2}
$$

Because the basis functions are orthonormal in the transformed representation, the overlap matrix in the Löwdin basis is the identity matrix $\mathbf{I}$, meaning that there are no shared off-diagonal overlap densities to partition. The gross electron population $N_\mu$ associated with the $\mu$-th Löwdin basis function is simply the diagonal element of the Löwdin density matrix

$$
N_\mu=P^{\text{Löwdin}}_{\mu\mu}=(\mathbf{S}^{1/2}\mathbf{P}\mathbf{S}^{1/2})_{\mu\mu}
$$

which is evaluated efficiently in our implementation by taking the dot product between the $\mu$-th row of the intermediate product $\mathbf{S}^{1/2}\mathbf{P}$ and the $\mu$-th row of the symmetric matrix $\mathbf{S}^{1/2}$. The total electron population $N_A$ assigned to atom $A$ is then obtained by summing the diagonal populations belonging to the basis functions on center $A$ as

$$
N_A=\sum_{\mu\in A}N_\mu
$$

using the mapping array `sys.bf2at`.

### 2. Net Atomic Charges

The net Löwdin atomic charge $q_A$ on atom $A$ is subsequently obtained by subtracting the assigned electron population from the nuclear charge according to

$$
q_A=Z_A-N_A
$$

which yields atomic partial charges that are significantly more robust against basis set enlargement than Mulliken charges.

---

## III. Mayer Bond Orders

While atomic partial charges provide information about the net distribution of electrons among individual atoms, they do not quantify the strength or covalent nature of chemical bonds connecting pairs of atoms. Mayer bond order analysis extends the concept of population analysis to interatomic pairs by evaluating the shared electron pair density between atoms directly within the non-orthogonal atomic orbital basis.

### 1. Mathematical Formulation

For a closed-shell electronic system described by the total density matrix $\mathbf{P}$ and the basis overlap matrix $\mathbf{S}$, the covalent bond order $B_{AB}$ between two distinct atoms $A$ and $B$ is defined by summing the products of the elements of the intermediate matrix $\mathbf{P}\mathbf{S}$ as

$$
B_{AB}=\sum_{\mu\in A}\sum_{\nu\in B}(\mathbf{P}\mathbf{S})_{\mu\nu}(\mathbf{P}\mathbf{S})_{\nu\mu}
$$

where the indices $\mu$ and $\nu$ run over all atomic basis functions centered on atoms $A$ and $B$, respectively. In our implementation, the matrix product $\mathbf{P}\mathbf{S}$ is evaluated first using dense matrix multiplication, and the pairwise off-diagonal products are accumulated for all atom pairs with $A\neq B$ using the basis-to-atom mapping array `sys.bf2at`.

### 2. Open-Shell and Generalized Systems

In spin-unrestricted or generalized electronic structure formalisms where electrons of different spins occupy distinct spatial orbitals, the total density is partitioned into alpha and beta spin components $\mathbf{P}^\alpha$ and $\mathbf{P}^\beta$. The Mayer bond order between atoms $A$ and $B$ is then calculated by evaluating the product matrices $\mathbf{P}^\alpha\mathbf{S}$ and $\mathbf{P}^\beta\mathbf{S}$ separately and accumulating their contributions according to

$$
B_{AB}=2\sum_{\mu\in A}\sum_{\nu\in B}\left[(\mathbf{P}^\alpha\mathbf{S})_{\mu\nu}(\mathbf{P}^\alpha\mathbf{S})_{\nu\mu}+(\mathbf{P}^\beta\mathbf{S})_{\mu\nu}(\mathbf{P}^\beta\mathbf{S})_{\nu\mu}\right]
$$

which incorporates spin polarization and yields bond orders that recover classical chemical valences for single, double, and triple bonds without requiring basis set orthogonalization.

---

## IV. Wiberg Bond Orders

The Wiberg bond index is an alternative measure of covalent bonding that was originally formulated in an orthogonal basis representation. Unlike the Mayer bond order, which operates directly on the non-orthogonal density and overlap matrices, the Wiberg bond order first projects the electronic density into the symmetrically orthogonalized Löwdin basis.

### 1. Symmetric Orthogonalization and Wiberg Index

To calculate Wiberg bond indices from non-orthogonal Gaussian basis functions, the overlap matrix $\mathbf{S}$ is diagonalized via eigenvalue decomposition $\mathbf{S}=\mathbf{U}\boldsymbol{\Lambda}\mathbf{U}^T$ to obtain the symmetric square root matrix $\mathbf{S}^{1/2}=\mathbf{U}\boldsymbol{\Lambda}^{1/2}\mathbf{U}^T$. The density matrix is then transformed into the orthonormal Löwdin basis as

$$
\mathbf{P}^{\text{ortho}}=\mathbf{S}^{1/2}\mathbf{P}\mathbf{S}^{1/2}
$$

which removes basis set overlap between distinct atomic centers. For a closed-shell system, the Wiberg bond order $W_{AB}$ between two distinct atoms $A$ and $B$ is calculated by summing the squares of the off-diagonal orthogonal density matrix elements as

$$
W_{AB}=\sum_{\mu\in A}\sum_{\nu\in B}(P^{\text{ortho}}_{\mu\nu})^2
$$

where the indices $\mu$ and $\nu$ run over the Löwdin basis functions mapped to atoms $A$ and $B$ via `sys.bf2at`.

### 2. Open-Shell and Generalized Systems

In generalized and spin-unrestricted calculations where the density matrix contains separate alpha and beta spin blocks $\mathbf{P}^\alpha$ and $\mathbf{P}^\beta$, each spin density matrix is transformed into the orthogonal Löwdin basis individually according to

$$
\mathbf{P}^{\alpha,\text{ortho}}=\mathbf{S}^{1/2}\mathbf{P}^\alpha\mathbf{S}^{1/2}
$$

and

$$
\mathbf{P}^{\beta,\text{ortho}}=\mathbf{S}^{1/2}\mathbf{P}^\beta\mathbf{S}^{1/2}
$$

and the total Wiberg bond index between atoms $A$ and $B$ is computed by summing the squared matrix elements across both spin manifolds as

$$
W_{AB}=2\sum_{\mu\in A}\sum_{\nu\in B}\left[(P^{\alpha,\text{ortho}}_{\mu\nu})^2+(P^{\beta,\text{ortho}}_{\mu\nu})^2\right]
$$

which provides a stable, basis-independent measure of electron sharing and chemical bond multiplicity between bonded atoms.
