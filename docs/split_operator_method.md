# Split-Operator Method

The Split-Operator (SPO) Fourier method is a highly efficient numerical algorithm developed by Feit, Fleck, and Steiger in 1982 to solve the time-dependent Schrödinger equation. It is widely used in quantum dynamics to propagate wavepackets on grid representations by dividing the time-evolution propagator into kinetic and potential energy components.

---

## I. Mathematical Formulation

### 1. Schrödinger Equation and Strang Splitting

The time evolution of a quantum state $|\psi(t)\rangle$ is governed by the time-dependent Schrödinger equation

$$
i\hbar\frac{d}{dt}|\psi(t)\rangle=\hat{H}|\psi(t)\rangle
$$

where $\hat{H}=\hat{T}+\hat{V}$ is the total Hamiltonian operator, consisting of the kinetic energy operator $\hat{T}$ and the potential energy operator $\hat{V}$. The formal solution for a small time step $\Delta t$ is

$$
|\psi(t+\Delta t)\rangle=\exp\left(-\frac{i}{\hbar}\hat{H}\Delta t\right)|\psi(t)\rangle
$$

Since the kinetic energy operator $\hat{T}$ and potential energy operator $\hat{V}$ do not commute, the exponential of their sum cannot be factored directly. The second-order Strang splitting scheme approximates the propagator as

$$
\exp\left(-\frac{i}{\hbar}\hat{H}\Delta t\right)=\exp\left(-\frac{i}{2\hbar}\hat{V}\Delta t\right)\exp\left(-\frac{i}{\hbar}\hat{T}\Delta t\right)\exp\left(-\frac{i}{2\hbar}\hat{V}\Delta t\right)+\mathcal{O}(\Delta t^3)
$$

which is accurate to second order in $\Delta t$.

---

## II. Representation and Multi-State Propagation

### 2. Position and Momentum Space Integrals

The split-operator method achieves high computational efficiency by evaluating each split propagator in the representation where that operator is diagonal. The potential energy operator $\hat{V}$ is diagonal in position space $\mathbf{r}$. Applying the potential propagator is a simple multiplication at each grid point as

$$
\psi'(\mathbf{r})=\exp\left(-\frac{i}{2\hbar}\mathbf{V}(\mathbf{r})\Delta t\right)\psi(\mathbf{r})
$$

In a multi-state nonadiabatic system, $\mathbf{V}(\mathbf{r})$ is an $N\times N$ matrix. The matrix exponential is evaluated by diagonalizing $\mathbf{V}(\mathbf{r})$ at each grid point using the unitary transformation matrix $\mathbf{U}(\mathbf{r})$ containing the eigenvectors and the diagonal matrix of adiabatic eigenvalues $\mathbf{W}(\mathbf{r})$ as

$$
\exp\left(-\frac{i}{2\hbar}\mathbf{V}(\mathbf{r})\Delta t\right)=\mathbf{U}(\mathbf{r})\exp\left(-\frac{i}{2\hbar}\mathbf{W}(\mathbf{r})\Delta t\right)\mathbf{U}^{\dagger}(\mathbf{r})
$$

which rotates the wavepacket into the adiabatic basis, applies the scalar phase factor updates, and rotates back. The kinetic energy operator $\hat{T}$ is diagonal in momentum space $\mathbf{p}$ (or wavenumber space $\mathbf{k}$). Applying the kinetic propagator involves transforming the wavefunction to momentum space using the Fast Fourier Transform (FFT), multiplying by the diagonal kinetic phase factors, and transforming back to position space using the Inverse Fast Fourier Transform (IFFT) as

$$
\tilde{\psi}''(\mathbf{k})=\exp\left(-\frac{i\hbar k^2}{2m}\Delta t\right)\tilde{\psi}'(\mathbf{k})
$$

where $m$ is the mass of the particle. The codebase implements these transforms by interfacing with the FFTW library, which generates optimized plans to carry out the multidimensional Fourier transforms.

---

## III. Boundary Conditions and Relaxation

### 3. Absorbing Boundary Potentials

To prevent unphysical reflections of the wavepacket at the grid boundaries, a complex absorbing potential (CAP) $-iV_{\text{cap}}(\mathbf{r})$ is added to the Hamiltonian. The effective Hamiltonian becomes

$$
\hat{H}_{\text{eff}}=\hat{T}+\hat{V}-iV_{\text{cap}}(\mathbf{r})
$$

which introduces a real exponential decay factor in the potential propagator that dampens the wavefunction as it approaches the boundaries of the grid, absorbing the outgoing flux.

### 4. Imaginary-Time Relaxation Trick

To find the ground-state wavefunction of a molecular system, the codebase implements the imaginary-time propagation trick. By substituting $\Delta t \to -i \Delta \tau$ where $\Delta \tau$ is a real parameter, the real oscillatory phase factors in the propagator turn into real decaying exponentials as

$$
\exp\left(-\frac{\hat{H}\Delta\tau}{\hbar}\right)
$$

This decay operator dampens high-energy eigenstates exponentially faster than the ground state. By repeatedly applying the split propagators and re-normalizing the wavepacket to unity at each step, the excited-state components vanish and the wavefunction relaxes to the exact numerical ground state.

---

## IV. Flux Analysis and Cross Section Calculations

### 5. Time-to-Energy Fourier Transform Flux Analysis

To compute energy-resolved reaction probabilities and scattering cross sections, the codebase implements a time-to-energy Fourier transform flux analysis. The energy-resolved wavefunction $\psi(E,\mathbf{r})$ is obtained from the time-propagated wavefunction $\psi(\mathbf{r},t)$ using the half-Fourier transform

$$
\psi(E,\mathbf{r})=\frac{1}{\sqrt{2\pi}}\int_0^{\infty}\psi(\mathbf{r},t)\exp\left(\frac{i}{\hbar}Et\right)dt
$$

which is discretized as the accumulated sum at each time step $\Delta t$

$$
A(E,\mathbf{r})=\sum_n\psi(\mathbf{r},t_n)\exp\left(\frac{i}{\hbar}Et_n\right)
$$

with $t_n=n\Delta t$. The energy-dependent reaction probability $P(E)$ is determined by integrating the quantum probability flux through a dividing surface normal to the scattering coordinate $d$ at position $x_s$ using the relation

$$
P(E)=\frac{\hbar}{\mu a_k(E)}\int\text{Im}\left[\psi^*(E,\mathbf{r})\nabla_d\psi(E,\mathbf{r})\right]d\mathbf{S}_d
$$

where $a_k(E)$ is the energy distribution of the initial wavepacket, $\mu$ is the effective mass of the incident coordinate, and the derivative $\nabla_d\psi(E,\mathbf{r})$ is evaluated spectrally in momentum space using Fast Fourier Transforms. In Cartesian coordinates, the scattering cross section $\sigma(E)$ is obtained from the reaction probability $P(E)$ by dividing by the transverse wavepacket density at the center of the coordinate system as

$$
\sigma(E)=\frac{P(E)}{r_{\text{perp}}}
$$

where the transverse normalization factor is defined as

$$
r_{\text{perp}}=\prod_{i=1}^{N-1}\sqrt{\frac{\gamma_i}{\pi}}
$$

with $\gamma_i$ representing the width parameters of the initial Gaussian wavepacket in the $N-1$ transverse directions.


### 6. Cylindrical Coordinate Transformation and Scaling

When modeling processes with cylindrical symmetry, such as diatomic collisions under the $J=0$ approximation, the radial coordinate $r$ introduces a $1/r$ term in the kinetic energy operator. Under this symmetry, the wavefunction has no dependence on the azimuthal angle $\theta$, reducing the 3D Schrödinger equation to a 2D problem in $(z, r)$. The Schrödinger equation is written as

$$
\hat{H}\Psi(r,z)=\left[-\frac{\hbar^2}{2\mu_z}\frac{\partial^2}{\partial z^2}-\frac{\hbar^2}{2\mu_r}\left(\frac{\partial^2}{\partial r^2}+\frac{1}{r}\frac{\partial}{\partial r}\right)+V(r,z)\right]\Psi(r,z)=E\Psi(r,z)
$$

where $\mu_z$ and $\mu_r$ are the coordinates' masses and $V(r,z)$ is the potential energy surface. To avoid non-Hermitian operators and allow the use of standard Cartesian Fast Fourier Transforms, the radial wavefunction is scaled using the relation

$$
\psi(r,z)=\sqrt{r}\Psi(r,z)
$$

where $r$ is the radial coordinate. We substitute this back into the Schrödinger equation by expressing the original wavefunction as $\Psi(r,z)=r^{-1/2}\psi(r,z)$. Differentiating $\Psi$ with respect to the radial coordinate $r$ gives the first derivative

$$
\frac{\partial\Psi}{\partial r}=-\frac{1}{2}r^{-3/2}\psi+r^{-1/2}\frac{\partial\psi}{\partial r}
$$

and the second derivative

$$
\frac{\partial^2\Psi}{\partial r^2}=\frac{3}{4}r^{-5/2}\psi-r^{-3/2}\frac{\partial\psi}{\partial r}+r^{-1/2}\frac{\partial^2\psi}{\partial r^2}
$$

which we substitute back into the radial kinetic energy operator term to yield the relation

$$
\left(\frac{\partial^2}{\partial r^2}+\frac{1}{r}\frac{\partial}{\partial r}\right)\Psi=r^{-1/2}\left(\frac{\partial^2\psi}{\partial r^2}+\frac{1}{4r^2}\psi\right)
$$

Multiplying the entire Schrödinger equation by $\sqrt{r}$ and substituting this relation simplifies the equation to a standard Cartesian form

$$
\left[-\frac{\hbar^2}{2\mu_z}\frac{\partial^2}{\partial z^2}-\frac{\hbar^2}{2\mu_r}\frac{\partial^2}{\partial r^2}+V_{\text{eff}}(r,z)\right]\psi(r,z)=E\psi(r,z)
$$

where the effective potential includes a centrifugal-like correction given by

$$
V_{\text{eff}}(r,z)=V(r,z)-\frac{\hbar^2}{8\mu_r r^2}
$$

with $\mu_r$ representing the mass associated with the radial coordinate. To satisfy the boundary condition at the origin where the wavefunction must vanish, the grid is extended symmetrically to negative radial values spanning $[-R_{\max}, R_{\max}]$, and the initial scaled wavefunction $\psi(r,z)$ is constructed with odd symmetry under $r\to-r$ as $\psi(r,z)=\text{sgn}(r)\sqrt{|r|}\Psi(r,z)$. This odd symmetry is preserved during propagation by the Cartesian kinetic energy operator, guaranteeing that the wavefunction remains zero at $r=0$. Because the negative coordinate region is a numerical extension, physical observables like position and momentum are computed using the absolute value of the radial grid coordinates. The scattering cross section $\sigma(E)$ is then computed by weighting the integrated flux with the cylindrical factor $\pi/\gamma_r$, where $\gamma_r$ is the radial wavepacket width parameter.

