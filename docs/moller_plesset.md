# Møller–Plesset Perturbation Theory

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of Møller–Plesset perturbation theory (MPPT) and how it is implemented in our scientific computing framework. If you want to understand how we can calculate the detailed correlation between electrons by treating it as a mathematical perturbation to the average Hartree–Fock field, this guide is written step-by-step for you.

---

## I. Perturbation Expansion and Spin Variants

To understand Møller–Plesset perturbation theory, we must first understand the limitations of the Hartree–Fock method. In Hartree–Fock, we assume that each electron moves in an average electric field created by all other electrons. In reality, electrons are negatively charged and repel each other individually. They perform a detailed dance to avoid getting close to one another, which is called electron correlation. Because Hartree–Fock misses this correlation, it underbinds molecules and gives inaccurate energies. Møller–Plesset perturbation theory is a way to calculate this correlation energy by starting with the Hartree–Fock solution and adding corrections to it step-by-step, treating the difference between the true electron interactions and the average field as a small mathematical perturbation.

### 1. Rayleigh–Schrödinger Perturbation Theory

To set up the perturbation theory, we divide the total electronic Hamiltonian operator $\hat{H}$ of the molecule into an unperturbed part $\hat{H}_0$, which is the sum of the one-electron Fock operators from Hartree–Fock, and a perturbation operator $\hat{V}$, which represents the difference between the true electron-electron repulsions and the average Hartree–Fock potential. We write this partition as

$$
\hat{H}=\hat{H}_0+\hat{V}
$$

where the reference state $| \Phi_0 \rangle$ is the Hartree–Fock Slater determinant. This reference determinant is an eigenstate of our unperturbed Hamiltonian $\hat{H}_0$, and its unperturbed energy is the sum of the orbital energies of the occupied electrons. For any excited Slater determinant $i$, the unperturbed energy is the sum of its occupied spin-orbital energies $\epsilon_p$ as

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

which describes how the perturbation mixes excited states into the wavefunction at higher orders. Our codebase implements this recurrence relation directly, enabling the user to calculate correlation energies up to any arbitrary perturbation order $k$ by recursively evaluating these coefficients and energy corrections.

### 2. Generalized and Spin-Unrestricted Møller–Plesset Variants

The specific equations of Møller–Plesset perturbation theory depend on the Hartree–Fock reference wavefunction we choose. Our codebase supports two main spin variants: Restricted Møller–Plesset (RMP) theory, which starts from a Restricted Hartree–Fock reference where all electrons are paired up in identical spatial orbitals; and Generalized Møller–Plesset (GMP) theory, which starts from a Generalized Hartree–Fock reference. GMP is formulated in a general spin-orbital basis where the spin components are allowed to mix, making it suitable for systems with non-collinear spin alignments or strong spin-orbit coupling.

---

## II. Nuclear Derivatives

### 3. Analytical Nuclear Gradient

To find the forces on the atoms in a molecule when using MP perturbation theory, we must calculate the derivative of the correlation energy with respect to the nuclear coordinates. This is more difficult than in Hartree–Fock because the Møller–Plesset wavefunction is not variationally optimized. This means that when the nuclei move, the molecular orbital coefficients change, and we must explicitly calculate their derivatives using the Coupled-Perturbed Hartree–Fock (CPHF) equations. By differentiating the energy expression of order $k$, the analytical nuclear gradient is given by

$$
\frac{dE^{(k)}}{dx}=\sum_{j\neq0}\left(\frac{dV_{0j}}{dx}C_j^{(k-1)}+V_{0j}\frac{dC_j^{(k-1)}}{dx}\right)
$$

where $x$ represents a nuclear coordinate, and the derivatives of the coefficients $dC_j^{(k-1)}/dx$ are calculated recursively. Evaluating this gradient requires propagating the derivatives of the molecular integrals and the molecular orbital coefficients, which is computationally demanding and complicated to implement.

### 4. Analytical Gradients via Dual Numbers

To simplify this process, our codebase implements an alternative approach using dual numbers. A dual number consists of a real part and a dual part containing a derivative, written as $a+b\epsilon$ where the dual unit satisfies $\epsilon^2=0$. Evaluating any smooth function $f$ on a dual number yields the function value and its exact derivative as

$$
f(x+dx\epsilon)=f(x)+f'(x)dx\epsilon
$$

which provides a way to calculate derivatives without numerical subtraction errors. To calculate the nuclear gradient of the correlation energy, we perturb the coordinates of the nuclei by adding the dual unit as $x\to x+\epsilon$ and compute the atomic orbital integrals as dual numbers containing their derivatives. By running the entire self-consistent field iterations, molecular orbital transformations, and perturbation calculations using the generic `ScalarDual(T)` type, the program automatically propagates the derivatives through the entire algorithm. The real part of the final result is the correlation energy, and the dual part is the exact analytical nuclear gradient. This dual number approach bypasses the need to write complex code to solve the CPHF response equations manually, reducing code complexity and preventing implementation errors.
