# Split-Operator Method

The Split-Operator (SPO) Fourier method is a highly efficient numerical algorithm developed by Feit, Fleck, and Steiger in 1982 to solve the time-dependent Schrödinger equation. It is widely used in quantum dynamics to propagate wavepackets on grid representations by dividing the time-evolution propagator into kinetic and potential energy components.

---

## 1. Mathematical Formulation

The time evolution of a quantum state $|\psi(t)\rangle$ is governed by the time-dependent Schrödinger equation

$$
i\hbar\frac{d}{dt}|\psi(t)\rangle=\hat{H}|\psi(t)\rangle
$$

where $\hat{H}=\hat{T}+\hat{V}$ is the total Hamiltonian operator, consisting of the kinetic energy operator $\hat{T}$ and the potential energy operator $\hat{V}$. The formal solution for a small time step $\Delta t$ is

$$
|\psi(t+\Delta t)\rangle=\exp\left(-\frac{i}{\hbar}\hat{H}\Delta t\right)|\psi(t)\rangle
$$

Since the kinetic energy operator $\hat{T}$ and potential energy operator $\hat{V}$ do not commute ($\left[\hat{T},\hat{V}\right]\neq0$), the exponential of their sum cannot be factored directly. The second-order Strang splitting scheme approximates the propagator as

$$
\exp\left(-\frac{i}{\hbar}\hat{H}\Delta t\right)=\exp\left(-\frac{i}{2\hbar}\hat{V}\Delta t\right)\exp\left(-\frac{i}{\hbar}\hat{T}\Delta t\right)\exp\left(-\frac{i}{2\hbar}\hat{V}\Delta t\right)+\mathcal{O}(\Delta t^3)
$$

which is accurate to second order in $\Delta t$.

---

## 2. Position and Momentum Space Representations

The split-operator method achieves high computational efficiency by evaluating each split propagator in the representation where that operator is diagonal.

The potential energy operator $\hat{V}$ is diagonal in position space $\mathbf{r}$. Applying the potential propagator is a simple multiplication at each grid point as

$$
\psi'(\mathbf{r})=\exp\left(-\frac{i}{2\hbar}\mathbf{V}(\mathbf{r})\Delta t\right)\psi(\mathbf{r})
$$

In a multi-state nonadiabatic system, $\mathbf{V}(\mathbf{r})$ is an $N\times N$ matrix. The matrix exponential is evaluated by diagonalizing $\mathbf{V}(\mathbf{r})$ at each grid point using the unitary transformation matrix $\mathbf{U}(\mathbf{r})$ containing the eigenvectors and the diagonal matrix of adiabatic eigenvalues $\mathbf{W}(\mathbf{r})$ as

$$
\exp\left(-\frac{i}{2\hbar}\mathbf{V}(\mathbf{r})\Delta t\right)=\mathbf{U}(\mathbf{r})\exp\left(-\frac{i}{2\hbar}\mathbf{W}(\mathbf{r})\Delta t\right)\mathbf{U}^{\dagger}(\mathbf{r})
$$

The kinetic energy operator $\hat{T}$ is diagonal in momentum space $\mathbf{p}$ (or wavenumber space $\mathbf{k}$). Applying the kinetic propagator involves transforming the wavefunction to momentum space using the Fast Fourier Transform (FFT), multiplying by the diagonal kinetic phase factors, and transforming back to position space using the Inverse Fast Fourier Transform (IFFT) as

$$
\tilde{\psi}''(\mathbf{k})=\exp\left(-\frac{i\hbar k^2}{2m}\Delta t\right)\tilde{\psi}'(\mathbf{k})
$$

where $m$ is the mass of the particle.

---

## 3. Absorbing Boundary Potentials

To prevent unphysical reflections of the wavepacket at the grid boundaries, a complex absorbing potential (CAP) $-iV_{\text{cap}}(\mathbf{r})$ is added to the Hamiltonian. The effective Hamiltonian becomes

$$
\hat{H}_{\text{eff}}=\hat{T}+\hat{V}-iV_{\text{cap}}(\mathbf{r})
$$

which introduces a real exponential decay factor in the potential propagator that dampens the wavefunction as it approaches the boundaries of the grid.
