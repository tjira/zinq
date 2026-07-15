# Ehrenfest Dynamics

Ehrenfest dynamics is a mixed quantum-classical molecular dynamics method that describes the interaction between classical nuclear degrees of freedom and quantum electronic degrees of freedom through a mean-field approximation. This approach partitions the system into a classical nuclear subsystem and a quantum electronic subsystem, allowing the simulation of nonadiabatic processes where multiple electronic states are involved.

---

## I. Quantum Electronic Propagation

### 1. Electronic Wavefunction Expansion

In the quantum description, the electronic wavefunction $\Psi(\mathbf{r},\mathbf{R}(t),t)$ is expanded in a diabatic basis $\{\phi_k(\mathbf{r})\}$ as

$$
\Psi(\mathbf{r},\mathbf{R}(t),t)=\sum_kc_k(t)\phi_k(\mathbf{r})
$$

where $c_k(t)$ are time-dependent complex coefficients and the basis functions do not depend on the nuclear coordinates, ensuring their derivatives with respect to the nuclear coordinates are zero. The electronic coefficients evolve according to the time-dependent Schrödinger equation

$$
i\hbar\frac{d}{dt}c_j(t)=\sum_kH_{jk}(\mathbf{R}(t))c_k(t)
$$

where $H_{jk}(\mathbf{R}(t))=\langle\phi_j|\hat{H}_{\text{el}}(\mathbf{R}(t))|\phi_k\rangle$ represents the electronic Hamiltonian matrix evaluated at the current classical nuclear coordinates. Under atomic units where $\hbar = 1$, the time derivative simplifies to

$$
\frac{dc_j(t)}{dt}=-i\sum_kH_{jk}(\mathbf{R}(t))c_k(t)
$$

which governs the propagation of the quantum subsystem along the classical path.

### 2. Implementation Tricks and Sub-stepping

Because electronic motion operates on a significantly faster timescale than classical nuclear motion, propagating both systems with the same classical time step $dt$ is numerically unstable. To circumvent this, the codebase implements a sub-stepping trick where the electronic coefficients are integrated over $N_{\text{steps}}$ sub-intervals of size $dt / N_{\text{steps}}$ during each classical step. The electronic TDSE is solved in the diabatic basis using the pre-allocated Runge–Kutta integrator (typically RK4), and the complex derivative is evaluated directly using optimized real-imaginary matrix-vector contractions to minimize execution overhead.

---

## II. Classical Nuclear Dynamics

### 3. Mean-Field Potential and Force

The classical nuclear degrees of freedom are governed by Newton's equations of motion

$$
M_A\frac{d^2\mathbf{R}_A}{dt^2}=\mathbf{F}_A
$$

where $M_A$ and $\mathbf{R}_A$ are the mass and position vector of nucleus $A$, and the classical force $\mathbf{F}_A$ is the negative gradient of the expectation value of the electronic Hamiltonian. Under the mean-field Ehrenfest approximation, the potential energy $E_{\text{pot}}$ is the expectation value of the electronic Hamiltonian

$$
E_{\text{pot}}=\langle\Psi|\hat{H}_{\text{el}}|\Psi\rangle
$$

which can be expressed in terms of the diabatic electronic coefficients as

$$
E_{\text{pot}}=\sum_{j,k}c_j^*(t)H_{jk}(\mathbf{R}(t))c_k(t)
$$

where the expectation value is guaranteed to be real because the Hamiltonian matrix is Hermitian. By defining the density matrix element $\rho_{jk}(t)=c_j^*(t)c_k(t)$, the mean-field potential energy becomes

$$
E_{\text{pot}}=\sum_{j,k}\text{Re}(c_j^*(t)c_k(t))H_{jk}(\mathbf{R}(t))
$$

which allows us to compute the mean-field force acting on nucleus $A$ as

$$
\mathbf{F}_A=-\nabla_{\mathbf{R}_A}E_{\text{pot}}=-\sum_{j,k}\text{Re}(c_j^*(t)c_k(t))\nabla_{\mathbf{R}_A}H_{jk}(\mathbf{R}(t))
$$

where the spatial gradients of the diabatic Hamiltonian matrix elements $\nabla_{\mathbf{R}_A} H_{jk}(\mathbf{R}(t))$ drive the classical nuclear trajectory.

---

## III. Basis Transformations and Populations

### 4. Initial State Setup and Projection

When the trajectory initial conditions are defined in the adiabatic basis, the initial wavefunction $|\Psi(0)\rangle$ is an eigenstate of the electronic Hamiltonian at the initial position $\mathbf{R}(0)$, which is written as

$$
\mathbf{H}(\mathbf{R}(0))\mathbf{u}_a=E_a\mathbf{u}_a
$$

where $\mathbf{u}_a$ is the eigenvector representing the active adiabatic state $a$. The initial electronic coefficients in the diabatic basis are set using the components of the eigenvector as

$$
c_j(0)=\mathbf{U}_{ja}
$$

where $\mathbf{U}$ is the unitary transformation matrix containing the eigenvectors of the Hamiltonian. If the simulation is performed in the adiabatic representation, the population fraction of the active adiabatic state $m$ at time $t$ is calculated by projecting the diabatic wavefunction back onto the adiabatic states as

$$
P_m(t)=\left|\sum_j\mathbf{U}_{jm}(\mathbf{R}(t))c_j(t)\right|^2
$$

where $\mathbf{U}(\mathbf{R}(t))$ is the unitary transformation matrix obtained by diagonalizing $\mathbf{H}(\mathbf{R}(t))$ at the current nuclear coordinates. When the simulation is performed in the diabatic representation, the population fraction of diabatic state $k$ is given by the squared magnitude of its coefficient

$$
P_k(t)=|c_k(t)|^2
$$

which completes the theoretical description of the electronic population dynamics.
