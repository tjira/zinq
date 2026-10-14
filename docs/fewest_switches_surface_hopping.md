# Fewest Switches Surface Hopping

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of Tully's Fewest Switches Surface Hopping (FSSH) method and how it is implemented in our framework. If you want to understand how we simulate quantum transitions by allowing classical trajectories to stochastically jump between different energy states, this guide is written step-by-step for you.

---

## I. Electronic Wavepacket Propagation

To understand Fewest Switches Surface Hopping, we must compare it to Ehrenfest dynamics. In Ehrenfest dynamics, the classical nuclei move on a single potential surface that is a weighted average of all electronic states. While this is simple, it can be unphysical because if a molecule splits into two products, one in the ground state and one in the excited state, Ehrenfest dynamics will predict a single product moving on a weird average state in between. Tully's surface hopping solves this by forcing each individual trajectory to move on a single active electronic state at a time. To simulate the quantum transitions, we propagate the electronic wavefunction along the trajectory and use its amplitudes to calculate a probability that the trajectory will suddenly jump (hop) to a different state.

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

To integrate the electronic coefficients, the program propagates the quantum state using tiny sub-steps during each classical step to ensure numerical stability. In the adiabatic representation, we obtain the eigenvectors of the Hamiltonian by diagonalizing the matrix at each step. However, diagonalization algorithms return eigenvectors with an arbitrary mathematical sign. If an eigenvector suddenly flips its sign from one step to the next, its overlap with the eigenvector from the previous step will become negative, causing the calculated derivative couplings to blow up. To prevent this, the code implements a phase-alignment check that computes the overlap of the eigenvectors between consecutive steps and flips the sign of the current eigenvector if the overlap is negative, keeping the phase smooth.

### 3. Time-Derivative Couplings and Approximations

In ab initio molecular dynamics using Configuration Interaction, the nonadiabatic time-derivative coupling matrix elements are evaluated directly from the nuclear velocity vector $\dot{\mathbf{R}}$ and the nonadiabatic derivative coupling vectors $\mathbf{d}_{ij}$ through the scalar product $\sigma_{ij}(t)=\dot{\mathbf{R}}\cdot\mathbf{d}_{ij}(t)$, where $\mathbf{d}_{ij}$ is computed from cross-geometry determinant overlaps as detailed in the Configuration–Interaction documentation. For model Hamiltonians, nonadiabatic coupling vectors $\mathbf{d}_{ij}$ are computed analytically in the coordinate-independent diabatic basis via the Hellmann–Feynman relation

$$
d_{ij,\alpha}=\frac{\sum_{k,l}U_{ki}U_{lj}\frac{\partial V_{kl}}{\partial R_\alpha}}{E_j-E_i}
$$

which directly yields the exact time-derivative coupling via $\sigma_{ij}(t)=\sum_\alpha\frac{P_\alpha}{M_\alpha}d_{ij,\alpha}$. Alternatively, when using wavefunction overlap coupling, the codebase approximates the time-derivative couplings using the finite-difference wavefunction overlap scheme proposed by Hammes-Schiffer and Tully as

$$
\sigma_{ij}(t)\approx\frac{S_{ji}-S_{ij}}{2\Delta t}
$$

where the overlap matrix elements are calculated from the unitary matrices containing the eigenvectors at consecutive steps as

$$
S_{ij}=\sum_kU_{ki}(t)U_{kj}(t-\Delta t)
$$

which propagates the electronic amplitudes using only the instantaneous eigenvalues and eigenvectors of the Hamiltonian.

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

When a trajectory hops from state $c$ to state $j$, the electronic potential energy changes by $\Delta E = E_j - E_c$. To conserve the total energy of the system, the classical kinetic energy of the nuclei must change by the exact same amount. In isotropic momentum adjustment, we achieve this by scaling the nuclear momentum vector $\mathbf{p}$ by a factor

$$
\gamma=\sqrt{\frac{E_{\text{kin}}-\Delta E}{E_{\text{kin}}}}
$$

where $E_{\text{kin}}$ is the initial kinetic energy of the nuclei. If the hop is uphill and the potential energy change is larger than the available kinetic energy, the term inside the square root is negative. In this case, the hop is rejected, and the trajectory remains on the active state $c$, which constitutes a frustrated hop. When nonadiabatic coupling vectors $\mathbf{d}_{cj}$ are explicitly available, momentum adjustment can equivalently be projected along the direction of $\mathbf{d}_{cj}$ to preserve angular momentum and align kinetic energy redistribution with the nonadiabatic transition coordinate.

---

## III. Decoherence Correction

### 6. Granucci–Persico Energy-Based Decoherence Correction

In standard Fewest Switches Surface Hopping, the electronic wavefunction retains phase coherence along the classical trajectory indefinitely, even after wavepackets propagating on distinct adiabatic surfaces have separated in nuclear phase space. This overcoherence can lead to spurious transitions and unphysical recrossings long after leaving the nonadiabatic interaction region. To restore physical wavepacket branching, our framework incorporates the energy-based decoherence correction proposed by Granucci and Persico, which is evaluated strictly in the adiabatic representation. For each inactive adiabatic state $j\ne c$, the characteristic decoherence time is defined as

$$
\tau_j=\frac{1+\alpha/E_{\text{kin}}}{|E_j-E_c|}
$$

where $\alpha$ is an empirical parameter, $E_{\text{kin}}$ is the nuclear kinetic energy, and $|E_j-E_c|$ is the adiabatic energy gap between state $j$ and the active state $c$. When $E_{\text{kin}}>10^{-4}\text{ a.u.}$, the electronic coefficients of all inactive states are damped at each integration substep $\Delta t$ according to

$$
c_j'(t)=c_j(t)\exp\left(-\frac{\Delta t}{2\tau_j}\right)
$$

while the active state coefficient is renormalized to conserve total electronic population through

$$
c_c'(t)=c_c(t)\sqrt{\frac{\max\left(0,1-\sum_{j\ne c}|c_j'(t)|^2\right)}{|c_c(t)|^2+10^{-14}}}
$$

which preserves the quantum phase of the active state while collapsing the electronic wavepacket toward the active adiabatic potential surface.
