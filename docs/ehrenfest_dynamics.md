# Ehrenfest Dynamics

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of Ehrenfest dynamics and how it is implemented in our framework. If you want to understand how we can simulate the coupled behavior of quantum electrons and classical nuclei using a mean-field approximation, this guide is written step-by-step for you.

---

## I. Quantum Electronic Propagation

To understand Ehrenfest dynamics, we must first look at the problem of simulating molecules that have been excited by light. In these systems, the electrons are not in their lowest-energy ground state, and they can transition between different electronic states. Because the nuclei are heavy and move slowly, we can treat them as classical particles. However, the electrons are light and fast, so we must treat them using quantum mechanics. Ehrenfest dynamics is a method that couples these two systems. The electronic state is described by a quantum wavefunction that changes over time, while the nuclei move classically. The force felt by the nuclei is a weighted average of the forces from all the electronic states, which is why it is called a mean-field method.

### 1. Electronic Wavefunction Expansion

In the quantum description, we expand the electronic wavefunction in a set of diabatic electronic states, which do not depend on the nuclear coordinates, as

$$
\Psi(\mathbf{r},\mathbf{R}(t),t)=\sum_kc_k(t)\phi_k(\mathbf{r})
$$

where $\mathbf{r}$ represents the electronic coordinates, $\mathbf{R}(t)$ represents the classical nuclear positions, $c_k(t)$ are time-dependent complex coefficients representing the amplitude of each electronic state, and $\phi_k(\mathbf{r})$ are the diabatic basis functions. The complex coefficients evolve in time according to the time-dependent Schrödinger equation

$$
i\hbar\frac{d}{dt}c_j(t)=\sum_kH_{jk}(\mathbf{R}(t))c_k(t)
$$

where $H_{jk}(\mathbf{R}(t))$ represents the elements of the electronic Hamiltonian matrix evaluated at the current classical nuclear coordinates. In atomic units where the reduced Planck constant is set to one, this equation simplifies to

$$
\frac{dc_j(t)}{dt}=-i\sum_kH_{jk}(\mathbf{R}(t))c_k(t)
$$

which describes how the electronic state vector changes as the nuclei move along their classical trajectory.

### 2. Implementation Tricks and Sub-stepping

Because electrons are extremely light, they move on a much faster timescale (femtoseconds or attoseconds) than the heavy nuclei (picoseconds). If we were to propagate both the electronic and nuclear equations using the same classical time step $dt$, the electronic calculation would become numerically unstable. To solve this, our codebase implements a sub-stepping procedure that divides the classical nuclear time step $dt$ into $N_{\text{steps}}$ smaller sub-intervals of size $dt / N_{\text{steps}}$. During each classical step, the program integrates the complex electronic coefficients over these small sub-steps using a fourth-order Runge–Kutta (RK4) integrator, while the nuclear coordinates are assumed to change linearly. This allows the quantum equations to remain stable without making the classical molecular dynamics propagation too expensive.

---

## II. Classical Nuclear Dynamics

### 3. Mean-Field Potential and Force

The movement of the classical nuclei is governed by Newton's second law of motion, which we write as

$$
M_A\frac{d^2\mathbf{R}_A}{dt^2}=\mathbf{F}_A
$$

where $M_A$ is the mass of nucleus $A$, $\mathbf{R}_A$ is its 3D position vector, and $\mathbf{F}_A$ is the classical force acting on it. In the diabatic representation, the potential energy $E_{\text{pot}}$ is the expectation value of the electronic Hamiltonian, which is calculated as

$$
E_{\text{pot}}=\sum_{j,k}c_j^*(t)H_{jk}(\mathbf{R}(t))c_k(t)
$$

where the asterisk denotes the complex conjugate. We can rewrite this potential energy in terms of the real part of the density matrix elements as

$$
E_{\text{pot}}=\sum_{j,k}\text{Re}(c_j^*(t)c_k(t))H_{jk}(\mathbf{R}(t))
$$

which is guaranteed to be a real number because the Hamiltonian matrix is Hermitian. The mean-field force acting on nucleus $A$ is the negative gradient of this potential energy, which is calculated as

$$
\mathbf{F}_A=-\sum_{j,k}\text{Re}(c_j^*(t)c_k(t))\nabla_{\mathbf{R}_A}H_{jk}(\mathbf{R}(t))
$$

where $\nabla_{\mathbf{R}_A}H_{jk}(\mathbf{R}(t))$ represents the spatial gradient of the Hamiltonian matrix elements, driving the classical nuclear motion on the average potential energy surface.

### 4. Adiabatic Formulation and Nonadiabatic Coupling

In ab initio molecular dynamics, electronic structure methods such as Configuration Interaction naturally yield adiabatic potential energy surfaces $E_k(\mathbf{R})$ and state gradients rather than diabatic states. In the adiabatic representation where the electronic state is expanded over adiabatic eigenstates with complex coefficients $a_k(t)$, the total electronic potential energy is the expectation value $E_{\text{pot}}=\sum_k |a_k(t)|^2 E_k(\mathbf{R})$. Applying the Hellmann–Feynman relation yields the mean-field adiabatic nuclear force acting on nucleus $A$ as

$$
\mathbf{F}_A=-\sum_k\rho_{kk}(t)\nabla_{\mathbf{R}_A}E_k(\mathbf{R})-\sum_{j\ne k}\text{Re}(\rho_{jk}(t))(E_k(\mathbf{R})-E_j(\mathbf{R}))\mathbf{d}_{jk,A}(\mathbf{R})
$$

where $\rho_{jk}(t)=a_j^*(t)a_k(t)$ is the electronic density matrix, $\nabla_{\mathbf{R}_A}E_k(\mathbf{R})$ is the adiabatic gradient of state $k$, and $\mathbf{d}_{jk,A}(\mathbf{R})=\langle\psi_j|\nabla_{\mathbf{R}_A}\psi_k\rangle$ is the nonadiabatic derivative coupling vector evaluated from cross-geometry determinant overlaps as detailed in the Configuration–Interaction documentation.

---

## III. Electronic State Representation and Populations

### 5. Initial State Setup and Populations

When conducting Ehrenfest dynamics directly in the adiabatic basis, the state amplitudes $a_k(t)$ represent the expansion coefficients in the local adiabatic eigenbasis. The initial conditions for an active adiabatic state $a$ are directly initialized as

$$
a_k(0)=\delta_{ka}
$$

and the population fraction of any adiabatic state $m$ at time $t$ is given by

$$
P_m(t)=|a_m(t)|^2
$$

without requiring coordinate-dependent basis transformations. When conducting dynamics in the diabatic basis, the electronic state is expanded in fixed diabatic basis functions with initial coefficients $c_j(0)=\delta_{ja}$, and the population fraction of diabatic state $k$ is given by

$$
P_k(t)=|c_k(t)|^2
$$

which describes how the electronic population shifts between the states over time.
