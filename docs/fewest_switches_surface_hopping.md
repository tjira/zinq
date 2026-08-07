# Fewest Switches Surface Hopping

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of Tully's Fewest Switches Surface Hopping (FSSH) method and how it is implemented in our framework. If you want to understand how we simulate quantum transitions by allowing classical trajectories to stochastically jump between different energy states, this guide is written step-by-step for you.

---

## I. Electronic Wavepacket Propagation

To understand Fewest Switches Surface Hopping, we must compare it to Ehrenfest dynamics. In Ehrenfest dynamics, the classical nuclei move on a single potential surface that is a weighted average of all electronic states. While this is simple, it can be unphysical: if a molecule splits into two products, one in the ground state and one in the excited state, Ehrenfest dynamics will predict a single product moving on a weird average state in between. Tully's surface hopping solves this by forcing each individual trajectory to move on a single active electronic state at a time. To simulate the quantum transitions, we propagate the electronic wavefunction along the trajectory and use its amplitudes to calculate a probability that the trajectory will suddenly jump (hop) to a different state.

### 1. Schrödinger Equation and Representation

We represent the electronic wavefunction as a linear combination of basis states as

$$
|\psi(t)\rangle=\sum_ic_i(t)|\phi_i\rangle
$$

where $c_i(t)$ are time-dependent complex coefficients representing the quantum amplitude of state $i$. In the adiabatic representation, the basis states are the eigenvectors of the electronic Hamiltonian at each nuclear coordinate. The coefficients in this basis evolve in time according to the adiabatic time-dependent Schrödinger equation

$$
\frac{dc_i(t)}{dt}=-iE_i(t)c_i(t)-\sum_j\sigma_{ij}(t)c_j(t)
$$

where $E_i(t)$ is the adiabatic energy of state $i$, and $\sigma_{ij}(t)=\dot{\mathbf{R}}\cdot\mathbf{d}_{ij}(t)$ is the nonadiabatic time-derivative coupling matrix element, which depends on the nuclear velocities $\dot{\mathbf{R}}$ and the nonadiabatic derivative coupling vectors $\mathbf{d}_{ij}(t)$. In the diabatic representation, the basis states do not change as the nuclei move, meaning their derivative couplings are zero. The coefficients in the diabatic basis evolve as

$$
\frac{dc_i(t)}{dt}=-i\sum_jH_{ij}(\mathbf{R}(t))c_j(t)
$$

where $H_{ij}(\mathbf{R}(t))$ represents the diabatic electronic Hamiltonian matrix elements at the current nuclear coordinates.

### 2. Implementation Tricks and Phase Alignment

To integrate the electronic coefficients, the program propagates the quantum state using tiny sub-steps during each classical step to ensure numerical stability. In the adiabatic representation, we obtain the eigenvectors of the Hamiltonian by diagonalizing the matrix at each step. However, diagonalization algorithms return eigenvectors with an arbitrary mathematical sign. If an eigenvector suddenly flips its sign from one step to the next, its overlap with the eigenvector from the previous step will become negative, causing the calculated derivative couplings to blow up. To prevent this, the code implements a phase-alignment check: it computes the overlap of the eigenvectors between consecutive steps and flips the sign of the current eigenvector if the overlap is negative, keeping the phase smooth.

### 3. Hammes-Schiffer–Tully Coupling Approximation

Calculating the analytical derivative coupling vectors $\mathbf{d}_{ij}$ is mathematically extremely difficult. To avoid this, our codebase approximates the time-derivative coupling elements $\sigma_{ij}$ using a finite-difference scheme proposed by Hammes-Schiffer and Tully, written as

$$
\sigma_{ij}(t)\approx\frac{S_{ji}-S_{ij}}{2\Delta t}
$$

where the overlap matrix elements are calculated from the unitary matrices containing the eigenvectors at consecutive steps as

$$
S_{ij}=\sum_kU_{ki}(t)U_{kj}(t-\Delta t)
$$

which allows us to propagate the wavepacket using only the eigenvalues and eigenvectors of the Hamiltonian.

---

## II. Tully Transition Probabilities and Hopping

### 4. Transition Probabilities

At each classical time step, the algorithm determines whether the trajectory should hop from the active state $c$ to a target state $j$. In the adiabatic basis, the probability of hopping from state $c$ to state $j$ during a time step $\Delta t$ is calculated using Tully's fewest switches formula as

$$
P_{c\to j}=\max\left(0,\frac{2\Delta t\text{Re}(c_c^{\ast}(t)c_j(t))\sigma_{cj}(t)}{\rho_{cc}(t)}\right)
$$

where $\rho_{cc}(t)=c_c(t)c_c^*(t)$ is the quantum population of the active state. In the diabatic basis, the hopping probability is calculated as

$$
P_{c\to j}=\max\left(0,\frac{2\Delta t\text{Im}(c_c(t)c_j^{\ast}(t)H_{cj}(\mathbf{R}(t)))}{\rho_{cc}(t)}\right)
$$

which drives the transitions. At each step, the program generates a random number between zero and one. If this random number is less than the calculated hopping probability, a transition is triggered.

### 5. Momentum Rescaling and Energy Conservation

When a trajectory hops from state $c$ to state $j$, the electronic potential energy changes by $\Delta E = E_j - E_c$. To conserve the total energy of the system, the classical kinetic energy of the nuclei must change by the exact same amount. We achieve this by scaling the nuclear momentum vector $\mathbf{p}$ by a factor

$$
\gamma=\sqrt{\frac{E_{\text{kin}}-\Delta E}{E_{\text{kin}}}}
$$

where $E_{\text{kin}}$ is the initial kinetic energy of the nuclei. If the hop is uphill and the potential energy change is larger than the available kinetic energy, the term inside the square root is negative. In this case, the hop is rejected, and the trajectory remains on the active state $c$, which constitutes a frustrated hop.
