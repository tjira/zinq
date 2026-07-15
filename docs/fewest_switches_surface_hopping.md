# Fewest Switches Surface Hopping

Fewest Switches Surface Hopping (FSSH) is a mixed quantum-classical molecular dynamics method developed by John Tully in 1990 to describe nonadiabatic processes in molecular systems. The method represents a trajectory-based approach where classical nuclei are propagated on a single active electronic state, while nonadiabatic transitions between electronic states are governed by stochastic hops determined by the evolution of the quantum electronic subsystem.

---

## I. Electronic Wavepacket Propagation

### 1. Schrödinger Equation and Representation

The electronic wavefunction $|\psi(t)\rangle$ is expanded in an orthonormal basis $\{|\phi_i\rangle\}$ as

$$
|\psi(t)\rangle=\sum_ic_i(t)|\phi_i\rangle
$$

where $c_i(t)$ are time-dependent complex coefficients. In the adiabatic representation, the basis states are eigenstates of the electronic Hamiltonian at each nuclear configuration, and the electronic coefficients evolve according to the adiabatic time-dependent Schrödinger equation

$$
\frac{dc_i(t)}{dt}=-iE_i(t)c_i(t)-\sum_j\sigma_{ij}(t)c_j(t)
$$

where $E_i(t)$ is the adiabatic energy of state $i$, and $\sigma_{ij}(t)=\langle\phi_i|\frac{\partial}{\partial t}|\phi_j\rangle=\dot{\mathbf{R}}\cdot\mathbf{d}_{ij}(t)$ is the nonadiabatic time-derivative coupling term. In the diabatic representation, the basis states do not depend on the nuclear coordinates, and the electronic coefficients evolve according to the diabatic time-dependent Schrödinger equation

$$
\frac{dc_i(t)}{dt}=-i\sum_jH_{ij}(\mathbf{R}(t))c_j(t)
$$

where $H_{ij}(\mathbf{R}(t))$ represents the diabatic electronic Hamiltonian matrix elements.

### 2. Implementation Tricks and Phase Alignment

To integrate the electronic coefficients, the codebase propagates the quantum subsystem using sub-steps of size $dt / N_{\text{steps}}$ during each classical step. In the adiabatic representation, standard diagonalization algorithms yield eigenvectors with arbitrary sign/phase. To prevent unphysical phase jumps between steps, the codebase implements a phase-alignment trick: it evaluates the overlap of corresponding eigenvectors at steps $t$ and $t-dt$, and flips the sign of eigenvector $j$ at step $t$ if the overlap is negative.

### 3. Hammes-Schiffer–Tully Coupling Approximation

Rather than requiring the explicit analytical calculation of non-adiabatic derivative coupling vectors $\mathbf{d}_{ij}$, the codebase approximates the time-derivative coupling matrix elements $\sigma_{ij}$ using the finite-difference Hammes-Schiffer–Tully scheme as

$$
\sigma_{ij}(t)\approx\frac{S_{ji}-S_{ij}}{2\Delta t}
$$

where the overlap terms are computed using the unitary eigenvectors from consecutive time steps as

$$
S_{ij}=\sum_kU_{ki}(t)U_{kj}(t-\Delta t)
$$

which permits efficient propagation using only energy eigenvalues and eigenvectors.

---

## II. Tully Transition Probabilities and Hopping

### 4. Transition Probabilities

At each classical time step, Tully's fewest switches algorithm determines whether a trajectory undergoes a transition from the active state $c$ to another state $j$. In the adiabatic basis, the transition probability is given by

$$
P_{c\to j}=\max\left(0,\frac{2\Delta t\text{Re}(c_c^{\ast}(t)c_j(t))\sigma_{cj}(t)}{\rho_{cc}(t)}\right)
$$

where $\rho_{cc}(t)=c_c(t)c_c^{\ast}(t)$ is the population of the active state. In the diabatic basis, the transition probability is given by

$$
P_{c\to j}=\max\left(0,\frac{2\Delta t\text{Im}(c_c(t)c_j^{\ast}(t)H_{cj}(\mathbf{R}(t)))}{\rho_{cc}(t)}\right)
$$

which drives the nonadiabatic transitions between diabatic states.

### 5. Momentum Rescaling and Energy Conservation

To conserve total energy during a nonadiabatic transition, the classical nuclear kinetic energy must compensate for the change in electronic potential energy. When a trajectory hops from state $c$ to state $j$, the change in potential energy is $\Delta E=E_j-E_c$. The nuclear momentum $\mathbf{p}$ is rescaled isotropically by a scaling factor

$$
\gamma=\sqrt{\frac{E_{\text{kin}}-\Delta E}{E_{\text{kin}}}}
$$

where $E_{\text{kin}}$ is the initial classical kinetic energy of the trajectory. If the transition is energetically uphill ($\Delta E>E_{\text{kin}}$), the final kinetic energy is negative, and the hop is rejected. The trajectory remains on the active state $c$, which constitutes a frustrated hop.
