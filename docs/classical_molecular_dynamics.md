# Classical Molecular Dynamics

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of classical molecular dynamics (MD) and how it is implemented in our framework. If you want to understand how we can simulate the physical movement of atoms in a molecule over time by solving Newton's equations of motion, this guide is written step-by-step for you.

---

## I. Trajectory Propagation

To understand molecular dynamics, we must look at how we partition a molecule. A molecule is made of light, fast-moving electrons and heavy, slow-moving atomic nuclei. Because the nuclei are thousands of times heavier than the electrons, we can make the simplifying assumption that the nuclei are treated as classical point-like particles that move according to Newton's laws of motion, while the electrons are treated quantum mechanically. The electrons adjust instantly to any movement of the nuclei, creating an average potential energy surface on which the nuclei move. Classical molecular dynamics is the process of calculating the forces acting on the nuclei at their current positions and moving them forward in time using small steps to simulate their physical trajectories.

### 1. Velocity Verlet Integration

To move the nuclei forward in time on a computer, we must integrate their equations of motion step-by-step. In our codebase, we use the Velocity Verlet algorithm, which is a second-order integration scheme. It is widely used in molecular simulations because it is symplectic, meaning it naturally conserves the total energy of the system over long times and does not suffer from numerical drift. For nuclear coordinates $\mathbf{r}$ and momenta $\mathbf{p}$ of a system with mass $M$, the coordinates are updated as

$$
\mathbf{r}(t+dt)=\mathbf{r}(t)+\frac{\mathbf{p}(t)}{M}dt+\frac{\mathbf{f}(t)}{2M}dt^2
$$

where $dt$ is the integration time step and $\mathbf{f}(t)$ represents the force vector acting on the nuclei at time $t$. The momenta are then updated in two separate half-steps. The first half-step update is calculated as

$$
\mathbf{p}\left(t+\frac{dt}{2}\right)=\mathbf{p}(t)+\frac{1}{2}\mathbf{f}(t)dt
$$

and the second half-step update, which is performed after calculating the new forces at the updated positions, is calculated as

$$
\mathbf{p}(t+dt)=\mathbf{p}\left(t+\frac{dt}{2}\right)+\frac{1}{2}\mathbf{f}(t+dt)dt
$$

which completes the step. This split update ensures that the coordinates and momenta are propagated in a highly stable, energy-conserving manner.

### 2. Force Evaluation on Potential Surfaces

The forces acting on the nuclei are computed as the negative spatial gradient of the potential energy surface. For a single-state adiabatic simulation where the molecule remains in its ground electronic state, the potential energy $V(\mathbf{r})$ is the electronic energy evaluated at the current geometry, yielding the classical force

$$
\mathbf{f}_i=-\frac{\partial V(\mathbf{r})}{\partial\mathbf{r}_i}
$$

for nuclear coordinates $\mathbf{r}_i$ of atom $i$, directing the nuclei towards regions of lower potential energy.

---

## II. Initial Condition Sampling

Generating realistic initial positions and momenta for molecular ensembles is essential for obtaining accurate thermodynamic averages and dynamical observables. In our codebase, initial conditions can be sampled according to classical thermal distributions or quantum phase space quasi-probability distributions.

### 3. Maxwell–Boltzmann Momentum Sampling

For classical thermal ensembles at a specified temperature $T$, nuclear coordinates are initially set to the reference equilibrium geometry $\mathbf{R}_0$, while nuclear momenta are sampled from the Maxwell–Boltzmann distribution

$$
P(\mathbf{P}_A)=\left(\frac{1}{2\pi M_A k_{\text{B}}T}\right)^{3/2}\exp\left(-\frac{\mathbf{P}_A^2}{2M_A k_{\text{B}}T}\right)
$$

where $M_A$ is the atomic mass of nucleus $A$ and $k_{\text{B}}$ is the Boltzmann constant. Each Cartesian momentum component is independently drawn from a univariate normal distribution with zero mean and variance $\sigma_{P,A}^2=M_A k_{\text{B}}T$. Following momentum sampling, the net linear momentum of the center of mass is subtracted to eliminate unphysical overall translational drift through space.

### 4. Wigner Distribution Sampling

In quantum molecular dynamics and nonadiabatic surface hopping, quantum nuclear zero-point energy and vibrational wavepacket delocalization are incorporated semiclassically by sampling initial positions and momenta from the Wigner distribution. For an ensemble of uncoupled quantum harmonic oscillators with normal mode frequencies $\omega_k$ at temperature $T$, the Wigner phase space density factorizes into independent Gaussian functions for each vibrational mode as

$$
\rho_{\text{W}}(Q_k,P_k)=\frac{1}{2\pi\sigma_{Q,k}\sigma_{P,k}}\exp\left(-\frac{Q_k^2}{2\sigma_{Q,k}^2}-\frac{P_k^2}{2\sigma_{P,k}^2}\right)
$$

where $Q_k$ and $P_k$ denote the mass-weighted normal coordinate and conjugate momentum. The coordinate variance is given by

$$
\sigma_{Q,k}^2=\frac{\hbar}{2\omega_k}\coth\left(\frac{\hbar\omega_k}{2k_{\text{B}}T}\right)
$$

and the conjugate momentum variance is given by

$$
\sigma_{P,k}^2=\frac{\hbar\omega_k}{2}\coth\left(\frac{\hbar\omega_k}{2k_{\text{B}}T}\right)
$$

which naturally recover the quantum harmonic zero-point fluctuations $\sigma_{Q,k}^2=\hbar/(2\omega_k)$ and $\sigma_{P,k}^2=\hbar\omega_k/2$ at zero temperature ($T=0$).

For analytical model potentials, Gaussian spatial width parameters $\boldsymbol{\gamma}$ can be directly provided, where $\gamma_k=\omega_k/\hbar$, defining variances $\sigma_{Q,k}^2=1/(2\gamma_k)$ and $\sigma_{P,k}^2=\hbar^2\gamma_k/2$.

For ab initio molecular systems, the vibrational frequencies and normal modes are determined from a harmonic vibrational frequency analysis. The mass-weighted Hessian matrix $\mathbf{H}^{\text{MW}}$ is evaluated at the optimized reference geometry $\mathbf{R}_0$ as

$$
H_{ia,jb}^{\text{MW}}=\frac{1}{\sqrt{M_a M_b}}\frac{\partial^2 V}{\partial R_{ia}\partial R_{jb}}
$$

and diagonalized to obtain normal mode frequencies $\omega_k=\sqrt{\lambda_k}$ and the orthonormal transformation matrix $\mathbf{L}$. Non-physical degrees of freedom corresponding to overall translation and rotation are discarded. Normal mode displacements $Q_k$ and momenta $P_k$ are independently sampled from Gaussian distributions with variances $\sigma_{Q,k}^2$ and $\sigma_{P,k}^2$, and transformed to Cartesian coordinates and momenta via

$$
\mathbf{R}=\mathbf{R}_0+\mathbf{M}^{-1/2}\mathbf{L}\mathbf{Q}
$$

and

$$
\mathbf{P}=\mathbf{M}^{1/2}\mathbf{L}\mathbf{P}_Q
$$

where $\mathbf{M}$ is the diagonal nuclear mass matrix. Center-of-mass linear and angular momenta are subsequently projected out to preserve rigid-body conservation laws during propagation.

---

## III. Non-Adiabatic Dynamics and Implementation

### 5. Ensemble Trajectory Management

In many simulations, running a single molecular trajectory is not enough to get reliable results because molecular behavior is statistical. To calculate observable properties (like average kinetic energy or chemical reaction rates), we run a large group of independent trajectories, which we call an ensemble. In our codebase, these trajectories are managed using an `Ensemble(T)` structure. This structure contains coordinate matrix `r`, momentum matrix `p`, and acceleration matrix `a` for all trajectories, along with a state vector `s` that stores the active electronic state index for each independent path.

### 6. Nonadiabatic Trajectories and Surface Hopping

When a molecular system transitions between distinct electronic potential energy surfaces, the single-surface Born–Oppenheimer approximation breaks down. Zinq supports nonadiabatic trajectory propagation using mean-field Ehrenfest dynamics as well as stochastic surface hopping. In Ehrenfest dynamics, the classical nuclei propagate on an effective potential energy surface obtained as an expectation value over the active electronic states, which is described in detail in the [Ehrenfest Dynamics](ehrenfest_dynamics.md) documentation. In surface hopping simulations, each trajectory in the ensemble evolves on a single active electronic state and undergoes stochastic hops mediated by transition probabilities, as detailed in the [Fewest Switches Surface Hopping](fewest_switches_surface_hopping.md) and [Landau–Zener Surface Hopping](landau_zener_surface_hopping.md) guides. When an electronic hop is accepted, the nuclear momentum vector $\mathbf{p}$ is rescaled to conserve total energy using the energy scaling factor

$$
\gamma=\sqrt{\frac{E_{\text{kin}}-\Delta E}{E_{\text{kin}}}}
$$

where $E_{\text{kin}}$ is the nuclear kinetic energy and $\Delta E = E_j - E_c$ is the change in electronic potential energy. If the hop is energetically uphill and the available kinetic energy is insufficient ($\Delta E > E_{\text{kin}}$), the transition is rejected as a frustrated hop and the trajectory continues on its active state.

### 7. Canonical Ensemble and Thermostats

To simulate molecular systems at constant temperature (the canonical NVT ensemble) rather than constant energy (the microcanonical NVE ensemble), our codebase implements temperature coupling through velocity rescaling and stochastic thermostats.

The Berendsen thermostat weakly couples the trajectory to an external thermal bath by scaling nuclear velocities at each step. The instantaneous temperature $T_{\text{inst}}$ is evaluated from the total nuclear kinetic energy $E_{\text{kin}}$ and the number of degrees of freedom $N_{\text{dof}}$ as

$$
T_{\text{inst}}=\frac{2E_{\text{kin}}}{N_{\text{dof}}}
$$

and the nuclear momenta are multiplied by a velocity rescaling factor

$$
\lambda=\sqrt{\max\left(0,1+\frac{\Delta t}{\tau}\left(\frac{T_{\text{target}}}{T_{\text{inst}}}-1\right)\right)}
$$

where $\tau$ is the characteristic coupling relaxation time and $T_{\text{target}}$ is the target bath temperature.

The Langevin thermostat incorporates thermal fluctuations and viscous damping by supplementing the physical forces with a friction force proportional to velocity and a stochastic Gaussian random force as

$$
M_A\frac{d^2\mathbf{R}_A}{dt^2}=\mathbf{f}_A-\gamma\mathbf{p}_A+\boldsymbol{\xi}_A(t)
$$

where $\gamma$ is the friction coefficient and $\boldsymbol{\xi}_A(t)$ is a stationary Gaussian white-noise vector satisfying the fluctuation-dissipation relation $\langle\xi_{Ai}(t)\xi_{Bj}(t')\rangle=2M_A\gamma T_{\text{target}}\delta_{AB}\delta_{ij}\delta(t-t')$, rigorously generating canonical ensemble sampling.
