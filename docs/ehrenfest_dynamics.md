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

Because electrons are extremely light, they move on a much faster timescale (femtoseconds or attoseconds) than the heavy nuclei (picoseconds). If we were to propagate both the electronic and nuclear equations using the same classical time step $dt$, the electronic calculation would become numerically unstable. To solve this, our codebase implements a sub-stepping procedure: we divide the classical nuclear time step $dt$ into $N_{\text{steps}}$ smaller sub-intervals of size $dt / N_{\text{steps}}$. During each classical step, the program integrates the complex electronic coefficients over these small sub-steps using a fourth-order Runge–Kutta (RK4) integrator, while the nuclear coordinates are assumed to change linearly. This allows the quantum equations to remain stable without making the classical molecular dynamics propagation too expensive.

---

## II. Classical Nuclear Dynamics

### 3. Mean-Field Potential and Force

The movement of the classical nuclei is governed by Newton's second law of motion, which we write as

$$
M_A\frac{d^2\mathbf{R}_A}{dt^2}=\mathbf{F}_A
$$

where $M_A$ is the mass of nucleus $A$, $\mathbf{R}_A$ is its 3D position vector, and $\mathbf{F}_A$ is the classical force acting on it. Under the Ehrenfest approximation, the potential energy $E_{\text{pot}}$ is the expectation value of the electronic Hamiltonian, which is calculated as

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

---

## III. Basis Transformations and Populations

### 4. Initial State Setup and Projection

At the start of the simulation, we often know the active state in the adiabatic representation, where the states are eigenstates of the electronic Hamiltonian at the initial nuclear position. To set up the initial conditions, we diagonalize the initial Hamiltonian matrix to solve the eigenvalue equation

$$
\mathbf{H}(\mathbf{R}(0))\mathbf{u}_a=E_a\mathbf{u}_a
$$

where $\mathbf{u}_a$ is the eigenvector representing the active adiabatic state $a$. The initial electronic coefficients in the diabatic basis are then set using the elements of the unitary transformation matrix $\mathbf{U}$ containing the eigenvectors as $c_j(0)=\mathbf{U}_{ja}$. During the dynamics, we can calculate the population fraction of any adiabatic state $m$ by projecting the diabatic coefficients back onto the adiabatic states as

$$
P_m(t)=\left|\sum_j\mathbf{U}_{jm}(\mathbf{R}(t))c_j(t)\right|^2
$$

where $\mathbf{U}(\mathbf{R}(t))$ is the unitary matrix of eigenvectors at the current coordinates. If we are running the simulation in the diabatic basis, the population fraction of diabatic state $k$ is simply the squared magnitude of its coefficient

$$
P_k(t)=|c_k(t)|^2
$$

which describes how the electronic population shifts between the states over time.
