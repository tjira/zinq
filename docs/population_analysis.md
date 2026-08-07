# Population Analysis

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of population analysis and how it is implemented in our scientific computing framework. If you want to understand how we can divide a continuous quantum mechanical electron cloud among the individual atoms in a molecule to calculate atomic charges, this guide is written step-by-step for you.

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
