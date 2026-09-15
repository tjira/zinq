# Møller–Plesset Perturbation Theory

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of Møller–Plesset perturbation theory (MPPT). If you want to understand how we can calculate the detailed correlation between electrons by treating it as a mathematical perturbation to the average Hartree–Fock field, this guide is written step-by-step for you.

---

## I. Perturbation Expansion

To understand Møller–Plesset perturbation theory, we must first understand the limitations of the Hartree–Fock method. In Hartree–Fock, we assume that each electron moves in an average electric field created by all other electrons. In reality, electrons are negatively charged and repel each other individually. They perform a detailed dance to avoid getting close to one another, which is called electron correlation. Because Hartree–Fock misses this correlation, it underbinds molecules and gives inaccurate energies. Møller–Plesset perturbation theory is a way to calculate this correlation energy by starting with the Hartree–Fock solution and adding corrections to it step-by-step, treating the difference between the true electron interactions and the average field as a small mathematical perturbation.

### 1. Rayleigh–Schrödinger Perturbation Theory

To set up the perturbation theory, we divide the total electronic Hamiltonian operator $\hat{H}$ of the molecule into an unperturbed part $\hat{H}_0$, which is the sum of the one-electron Fock operators from Hartree–Fock, and a perturbation operator $\hat{V}$, which represents the difference between the true electron-electron repulsions and the average Hartree–Fock potential. We write this partition as

$$
\hat{H}=\hat{H}_0+\hat{V}
$$

where the reference state $|\Phi_0\rangle$ is the Hartree–Fock Slater determinant. This reference determinant is an eigenstate of our unperturbed Hamiltonian $\hat{H}_0$, and its unperturbed energy is the sum of the orbital energies of the occupied electrons. For any excited Slater determinant $i$, the unperturbed energy is the sum of its occupied spin-orbital energies $\epsilon_p$ as

$$
E_i^{(0)}=\sum_{p\in\text{det}_i}\epsilon_p
$$

where the sum runs over the indices of the spin-orbitals occupied in determinant $i$. Rayleigh–Schrödinger perturbation theory assumes that we can expand the exact energy and wavefunction of the system as power series in terms of the perturbation. Using the intermediate normalization convention where the overlap of the unperturbed ground state with the exact wavefunction is exactly one, the electronic energy correction at order $k$ is given by

$$
E^{(k)}=\sum_{j\neq0}V_{0j}C_j^{(k-1)}
$$

for any order $k\ge2$, where the first-order energy correction is the expectation value $E^{(1)}=V_{00}$, and $V_{ij}=\langle\Phi_i|\hat{V}|\Phi_j\rangle$ represents the perturbation matrix elements between Slater determinants $i$ and $j$. The corresponding first-order wavefunction coefficients $C_i^{(1)}$ for any excited determinant $i\neq0$ are given by

$$
C_i^{(1)}=\frac{V_{i0}}{E_0^{(0)}-E_i^{(0)}}
$$

where the denominator $E_0^{(0)}-E_i^{(0)}$ represents the excitation energy of determinant $i$ from the ground state. For higher orders $k\ge2$, the coefficients are updated recursively using the formula

$$
C_i^{(k)}=\frac{\sum_{j\neq0}V_{ij}C_j^{(k-1)}-\sum_{j=1}^{k-1}E^{(j)}C_i^{(k-j)}}{E_0^{(0)}-E_i^{(0)}}
$$

which describes how the perturbation mixes excited states into the wavefunction at higher orders. This recurrence relation enables the calculation of correlation energies up to any arbitrary perturbation order $k$ by recursively evaluating these coefficients and energy corrections.

---

## II. Nuclear Derivatives

### 2. Analytical Nuclear Gradient

To find the forces on the atoms in a molecule when using MP perturbation theory, we must calculate the derivative of the total energy with respect to the nuclear coordinates. Because the Møller–Plesset wavefunction is not variationally optimized, its energy is not stationary with respect to changes in the molecular orbital coefficients, meaning the nuclear gradient depends explicitly on the first-order response of the molecular orbitals and their energies. The total analytical nuclear gradient at perturbation order $k$ combines the ground-state Hartree–Fock gradient and the correlation energy gradients as

$$
\frac{dE_{\text{tot}}}{dx}=\frac{dE_{\text{HF}}}{dx}+\sum_{m=2}^k\frac{dE^{(m)}}{dx}
$$

where the derivative of each correlation correction $E^{(m)}$ is obtained by differentiating the Rayleigh–Schrödinger energy expression as

$$
\frac{dE^{(m)}}{dx}=\sum_{j\neq0}\left(\frac{dV_{0j}}{dx}C_j^{(m-1)}+V_{0j}\frac{dC_j^{(m-1)}}{dx}\right)
$$

which involves the derivatives of the perturbation matrix elements $\frac{dV_{0j}}{dx}$ and the recursive differentiation of the wavefunction coefficients $\frac{dC_j^{(m-1)}}{dx}$ through the orbital energy denominators.

### 3. Differentiation via CPHF and Dual Numbers

Evaluating the correlation gradient directly by differentiating the determinant matrix elements and recursive perturbation equations by hand is mathematically tedious and prone to algebraic errors. To resolve this, mean-field response theory is combined with forward-mode automatic differentiation using dual numbers.

First, the Coupled-Perturbed Hartree–Fock (CPHF) equations are solved to determine the first-order mean-field response to nuclear displacements. Specifically, CPHF yields the molecular orbital coefficient derivatives $\frac{d\mathbf{C}}{dx}$ and the orbital energy derivatives $\frac{d\boldsymbol{\epsilon}}{dx}$, which describe how the molecular orbitals rotate and shift in energy as the atoms move.

Second, rather than manually differentiating the four-index molecular orbital transformation, the Slater–Condon rules, and the recursive Rayleigh–Schrödinger equations, the post-Hartree–Fock perturbation series is differentiated using forward-mode automatic differentiation. Seeding the molecular orbital coefficients, orbital energies, and two-electron repulsion integrals with their exact first-order derivatives ($\frac{d\mathbf{C}}{dx}$, $\frac{d\boldsymbol{\epsilon}}{dx}$, and $\frac{d\mathbf{g}}{dx}$) automatically propagates the exact chain-rule derivatives through the perturbation series using dual numbers. Evaluating the perturbation series yields the correlation energy $E^{(m)}$ in the real component and the exact analytical correlation gradient $\frac{dE^{(m)}}{dx}$ in the dual component, completing the total analytical nuclear gradient.
