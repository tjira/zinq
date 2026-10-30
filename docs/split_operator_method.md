# Split-Operator Method

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of the Split-Operator (SPO) Fourier method and how it is implemented in our quantum dynamics codebase. If you want to understand how we simulate the movements of a quantum wavepacket in time by switching back and forth between position and momentum space using Fourier transforms, this guide is written step-by-step for you.

---

## I. Mathematical Formulation

To understand the Split-Operator method, we must first look at what a quantum wavepacket is. In quantum mechanics, a particle (such as an atom or a molecule) does not have a single exact position. Instead, it is described by a wavefunction, which is a complex-valued wave that represents the probability distribution of where the particle might be found. The way this wavefunction changes over time is governed by the time-dependent Schrödinger equation. The main difficulty in solving this equation is that the energy operator, called the Hamiltonian, contains two parts, namely the potential energy (which depends on the positions of the atoms) and the kinetic energy (which depends on their velocities or momenta). Because these two operators do not commute, we cannot calculate them easily in the same representation. The Split-Operator method solves this by splitting the time step into small parts, allowing us to calculate each energy component in the representation where it is easiest to evaluate.

### 1. Schrödinger Equation and Strang Splitting

The time evolution of a quantum state wavefunction $|\psi(t)\rangle$ is governed by the time-dependent Schrödinger equation

$$
i\hbar\frac{d}{dt}|\psi(t)\rangle=\hat{H}|\psi(t)\rangle
$$

where $\hat{H}=\hat{T}+\hat{V}$ is the Hamiltonian operator, consisting of the kinetic energy operator $\hat{T}$ and the potential energy operator $\hat{V}$. The formal solution for how the wavefunction changes over a small time step $\Delta t$ is written as

$$
|\psi(t+\Delta t)\rangle=\exp\left(-\frac{i}{\hbar}\hat{H}\Delta t\right)|\psi(t)\rangle
$$

where the term containing the exponential is the propagator. Because the kinetic energy operator $\hat{T}$ (which involves derivatives) and the potential energy operator $\hat{V}$ (which involves coordinate values) do not commute, the exponential of their sum cannot be factored directly. To solve this, we use the second-order Strang splitting scheme to approximate the propagator as

$$
\exp\left(-\frac{i}{\hbar}\hat{H}\Delta t\right)=\exp\left(-\frac{i}{2\hbar}\hat{V}\Delta t\right)\exp\left(-\frac{i}{\hbar}\hat{T}\Delta t\right)\exp\left(-\frac{i}{2\hbar}\hat{V}\Delta t\right)+\mathcal{O}(\Delta t^3)
$$

where we apply half of the potential energy propagation, then the full kinetic energy propagation, and finally the second half of the potential energy propagation. This symmetric splitting cancels out the first-order error terms, making the integration highly accurate.

---

## II. Representation and Multi-State Propagation

### 2. Position and Momentum Space Integrals

The main advantage of the Split-Operator method is that it evaluates each split propagator in the representation where that operator is diagonal. The potential energy operator $\hat{V}$ is diagonal in position space. In this space, the potential propagator is just a simple multiplication at each grid point, written as

$$
\psi'(\mathbf{r})=\exp\left(-\frac{i}{2\hbar}\mathbf{V}(\mathbf{r})\Delta t\right)\psi(\mathbf{r})
$$

where $\mathbf{r}$ represents the position coordinates on our grid. If the system has multiple electronic states, the potential energy $\mathbf{V}(\mathbf{r})$ at each grid point is a matrix rather than a single number. To apply the propagator to this matrix, the program diagonalizes the potential matrix at each grid point using the unitary transformation matrix $\mathbf{U}(\mathbf{r})$ containing its eigenvectors and the diagonal matrix of adiabatic eigenvalues $\mathbf{W}(\mathbf{r})$ as

$$
\exp\left(-\frac{i}{2\hbar}\mathbf{V}(\mathbf{r})\Delta t\right)=\mathbf{U}(\mathbf{r})\exp\left(-\frac{i}{2\hbar}\mathbf{W}(\mathbf{r})\Delta t\right)\mathbf{U}^{\dagger}(\mathbf{r})
$$

which rotates the wavefunction into the adiabatic basis, applies the simple phase updates, and rotates it back to the original diabatic basis. The kinetic energy operator $\hat{T}$ is diagonal in momentum space $\mathbf{p}$. To apply the kinetic propagator, the program uses the Fast Fourier Transform (FFT) to convert the wavefunction from position space to momentum space, multiplies the values by the kinetic phase factors, and then uses the Inverse Fast Fourier Transform (IFFT) to convert the wavefunction back to position space. The momentum space update is calculated as

$$
\tilde{\psi}''(\mathbf{k})=\exp\left(-\frac{i\hbar k^2}{2m}\Delta t\right)\tilde{\psi}'(\mathbf{k})
$$

where $m$ is the mass of the particle and $\mathbf{k}$ is the wavevector representing momentum. Our codebase implements these transforms using a native multi-dimensional radix-2 Cooley-Tukey fast Fourier transform algorithm with bit-reversal indexing and stride-based dimensional decomposition.

---

## III. Boundary Conditions and Relaxation

### 3. Absorbing Boundary Potentials

Because we must run our simulations on a finite grid of points, we face a physical problem where the wavepacket reaches the edge of the grid and reflects off the boundary, traveling backward and interfering with itself. This is unphysical, as a real wavepacket would simply fly off into space. To prevent these reflections, we add an imaginary absorbing potential $-iV_{\text{cap}}(\mathbf{r})$ near the edges of the grid. The effective Hamiltonian becomes

$$
\hat{H}_{\text{eff}}=\hat{T}+\hat{V}-iV_{\text{cap}}(\mathbf{r})
$$

where the imaginary term acts like a sponge, introducing a real decaying exponential factor in the potential propagator that dampens the wavefunction to zero as it approaches the boundaries of the grid, absorbing the outgoing flux.

### 4. Imaginary-Time Relaxation Trick

If we want to find the lowest-energy ground-state wavefunction of a system, we can use a mathematical trick called imaginary-time propagation. By substituting the real time step with an imaginary time step $\Delta t \to -i \Delta \tau$, the oscillatory phase factors in our propagator turn into real decaying exponentials as

$$
\exp\left(-\frac{\hat{H}\Delta\tau}{\hbar}\right)
$$

which decays the amplitudes of the different energy states. Because higher-energy excited states decay exponentially faster than the lowest-energy ground state, repeatedly applying this propagator and re-normalizing the total probability of the wavefunction back to one will cause all the excited states to disappear, leaving only the exact numerical ground state.

---

## IV. Flux Analysis and Cross Section Calculations

### 5. Time-to-Energy Fourier Transform Flux Analysis

To calculate chemical reaction rates and scattering cross sections, we want to know how much of the wavepacket passes through a specific plane (called a dividing surface) at different energies. We first calculate the energy-resolved wavefunction $\psi(E,\mathbf{r})$ by taking the Fourier transform of the time-propagated wavefunction $\psi(\mathbf{r},t)$ as

$$
\psi(E,\mathbf{r})=\frac{1}{\sqrt{2\pi}}\int_0^{\infty}\psi(\mathbf{r},t)\exp\left(\frac{i}{\hbar}Et\right)dt
$$

which we discretize on our computer as an accumulated sum over all time steps $t_n=n\Delta t$ as

$$
A(E,\mathbf{r})=\sum_n\psi(\mathbf{r},t_n)\exp\left(\frac{i}{\hbar}Et_n\right)
$$

where the sum accumulates the wavepacket values. The reaction probability $P(E)$ is then computed by integrating the quantum probability flux through the dividing surface as

$$
P(E)=\frac{\hbar}{\mu a_k(E)}\int\text{Im}\left[\psi^*(E,\mathbf{r})\nabla_d\psi(E,\mathbf{r})\right]d\mathbf{S}_d
$$

where $a_k(E)$ is the energy distribution of the initial wavepacket, $\mu$ is the mass, and the spatial derivative $\nabla_d\psi(E,\mathbf{r})$ along the coordinate normal to the surface is calculated in momentum space using FFTs. The scattering cross section $\sigma(E)$ is then calculated by dividing the reaction probability by the transverse density of the wavepacket as

$$
\sigma(E)=\frac{P(E)}{r_{\text{perp}}}
$$

where the normalization factor is defined as

$$
r_{\text{perp}}=\prod_{i=1}^{N-1}\sqrt{\frac{\gamma_i}{\pi}}
$$

with $\gamma_i$ representing the width parameters of the initial wavepacket.

### 6. Cylindrical Coordinate Transformation and Scaling

When a collision has cylindrical symmetry (like a diatomic molecule colliding with an atom), we can reduce the 3D Schrödinger equation to a 2D problem in coordinates $(z, r)$, where $z$ is the axis and $r$ is the radial distance. The Hamiltonian in these coordinates is written as

$$
\hat{H}\Psi(r,z)=\left[-\frac{\hbar^2}{2\mu_z}\frac{\partial^2}{\partial z^2}-\frac{\hbar^2}{2\mu_r}\left(\frac{\partial^2}{\partial r^2}+\frac{1}{r}\frac{\partial}{\partial r}\right)+V(r,z)\right]\Psi(r,z)=E\Psi(r,z)
$$

where the radial term contains a first derivative $1/r$ that is not Hermitian and cannot be solved using standard Cartesian FFTs. To solve this, we define a scaled wavefunction as

$$
\psi(r,z)=\sqrt{r}\Psi(r,z)
$$

and substitute this back into the equations. Differentiating this scaled wavefunction gives the first radial derivative as

$$
\frac{\partial\Psi}{\partial r}=-\frac{1}{2}r^{-3/2}\psi+r^{-1/2}\frac{\partial\psi}{\partial r}
$$

and the second radial derivative as

$$
\frac{\partial^2\Psi}{\partial r^2}=\frac{3}{4}r^{-5/2}\psi-r^{-3/2}\frac{\partial\psi}{\partial r}+r^{-1/2}\frac{\partial^2\psi}{\partial r^2}
$$

which we combine to simplify the radial kinetic energy operator to

$$
\left(\frac{\partial^2}{\partial r^2}+\frac{1}{r}\frac{\partial}{\partial r}\right)\Psi=r^{-1/2}\left(\frac{\partial^2\psi}{\partial r^2}+\frac{1}{4r^2}\psi\right)
$$

where the first-derivative terms cancel out. Multiplying the entire Schrödinger equation by $\sqrt{r}$ yields a standard Cartesian-like equation

$$
\left[-\frac{\hbar^2}{2\mu_z}\frac{\partial^2}{\partial z^2}-\frac{\hbar^2}{2\mu_r}\frac{\partial^2}{\partial r^2}+V_{\text{eff}}(r,z)\right]\psi(r,z)=E\psi(r,z)
$$

where the effective potential includes a centrifugal correction term calculated as

$$
V_{\text{eff}}(r,z)=V(r,z)-\frac{\hbar^2}{8\mu_r r^2}
$$

which can be solved using standard FFT algorithms. To satisfy the physical boundary condition where the wavefunction must be zero at $r=0$, the radial grid is extended symmetrically to negative values spanning $[-R_{\max}, R_{\max}]$, and the initial wavefunction is constructed with odd symmetry as $\psi(r,z)=\text{sgn}(r)\sqrt{|r|}\Psi(r,z)$. The Cartesian kinetic energy propagator preserves this odd symmetry, guaranteeing that the wavefunction remains zero at $r=0$.

### 7. Partial Wave Expansion and Total Reaction Cross Section

To determine total integral reactive scattering cross sections in three-dimensional collisional dynamics, the scattering wavefunction is decomposed into partial waves indexed by the total angular momentum quantum number $J$. For each partial wave channel $J$, the effective potential on the radial coordinate $r$ receives an additional orbital centrifugal barrier defined as

$$
V_{\text{cent}}(r)=\frac{\hbar^2 J(J+1)}{2\mu r^2}
$$

where $\mu$ is the collision reduced mass along the radial reaction coordinate and $\hbar=1$ in atomic units. This centrifugal barrier suppresses reactivity at large impact parameters corresponding to high angular momenta, ensuring that the partial wave series converges within a finite number of partial waves $J \le J_{\max}$.

The maximum angular momentum $J_{\max}$ required to achieve convergence can either be supplied explicitly or estimated automatically from the maximum collision energy $E_{\max}$ and the potential barrier peak height $V_{\text{peak}}$ located at $r_{\text{peak}}$ as

$$
J_{\max}=\left\lceil\sqrt{2\mu\max(1.1E_{\max}-V_{\text{peak}}, 0)}\,r_{\text{peak}}\right\rceil+4
$$

which guarantees that all open classical impact parameters capable of overcoming the effective potential barrier are captured in the dynamics.

For each angular momentum $J$, the time-dependent wavepacket is propagated on the spatial grid using the split-operator method with the centrifugal term added directly to the diagonal potential energy operator. The energy-resolved reaction probability $P_J(E)$ is extracted at the dividing surface via time-to-energy Fourier transform flux analysis. The total integral reaction cross section $\sigma(E)$ is then computed by summing the degeneracy-weighted partial wave contributions as

$$
\sigma(E)=\frac{\pi}{2\mu E}\sum_{J=0}^{J_{\max}}(2J+1)\Delta J\,P_J(E)
$$

where $2J+1$ is the spatial degeneracy factor associated with the magnetic quantum number $M_J$, $E$ is the asymptotic collision energy measured relative to the reactant threshold, and $\Delta J$ is the angular momentum sampling step. Because each partial wave simulation is completely decoupled from the others, the evaluation over angular momentum channels is parallelized dynamically across multiple worker threads using atomic counters, yielding high computational throughput while preserving full quantum coherence within each partial wave channel.
