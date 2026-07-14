# Fewest Switches Surface Hopping

Fewest Switches Surface Hopping (FSSH) is a mixed quantum-classical molecular dynamics method developed by John Tully in 1990 to describe nonadiabatic processes in molecular systems. The method represents a trajectory-based approach where classical nuclei are propagated on a single active electronic state, while nonadiabatic transitions between electronic states are governed by stochastic hops determined by the evolution of the quantum electronic subsystem.

---

## 1. Electronic Wavepacket Propagation

The electronic wavefunction $|\psi(t)\rangle$ is expanded in an orthonormal basis $\{|\phi_i\rangle\}$ as

$$
|\psi(t)\rangle=\sum_ic_i(t)|\phi_i\rangle
$$

where $c_i(t)$ are time-dependent complex coefficients. 

In the adiabatic representation, the basis states are eigenstates of the electronic Hamiltonian at each nuclear configuration, and the electronic coefficients evolve according to the adiabatic time-dependent Schrödinger equation

$$
\frac{dc_i(t)}{dt}=-iE_i(t)c_i(t)-\sum_j\sigma_{ij}(t)c_j(t)
$$

where $E_i(t)$ is the adiabatic energy of state $i$, and $\sigma_{ij}(t)=\langle\phi_i|\frac{\partial}{\partial t}|\phi_j\rangle=\dot{\mathbf{R}}\cdot\mathbf{d}_{ij}(t)$ is the nonadiabatic time-derivative coupling term. The nonadiabatic coupling terms are computed using the Hammes-Schiffer–Tully approximation as

$$
\sigma_{ij}(t)\approx\frac{\langle\phi_j(t)|\phi_i(t-\Delta t)\rangle-\langle\phi_i(t)|\phi_j(t-\Delta t)\rangle}{2\Delta t}
$$

which is evaluated via the overlap of the electronic transformation eigenvectors at subsequent time steps.

In the diabatic representation, the basis states do not depend on the nuclear coordinates, and the electronic coefficients evolve according to the diabatic time-dependent Schrödinger equation

$$
\frac{dc_i(t)}{dt}=-i\sum_jH_{ij}(\mathbf{R}(t))c_j(t)
$$

where $H_{ij}(\mathbf{R}(t))$ represents the diabatic electronic Hamiltonian matrix elements.

---

## 2. Tully Transition Probabilities

At each classical time step, Tully's fewest switches algorithm determines whether a trajectory undergoes a transition from the active state $c$ to another state $j$.

In the adiabatic basis, the transition probability is given by

$$
P_{c\to j}=\max\left(0,\frac{2\Delta t\text{Re}(c_c^{\ast}(t)c_j(t))\sigma_{cj}(t)}{\rho_{cc}(t)}\right)
$$

where $\rho_{cc}(t)=c_c(t)c_c^{\ast}(t)$ is the population of the active state.

In the diabatic basis, the transition probability is given by

$$
P_{c\to j}=\max\left(0,\frac{2\Delta t\text{Im}(c_c(t)c_j^{\ast}(t)H_{cj}(\mathbf{R}(t)))}{\rho_{cc}(t)}\right)
$$

which drives the nonadiabatic transitions between diabatic states.

---

## 3. Momentum Rescaling and Energy Conservation

To conserve total energy during a nonadiabatic transition, the classical nuclear kinetic energy must compensate for the change in electronic potential energy. When a trajectory hops from state $c$ to state $j$, the change in potential energy is

$$
\Delta E=E_j-E_c
$$

where $E_c$ and $E_j$ are the potential energies of the initial and final states. The nuclear momentum $\mathbf{P}$ is rescaled isotropically by a scaling factor

$$
\gamma=\sqrt{\frac{E_{\text{kin}}-\Delta E}{E_{\text{kin}}}}
$$

where $E_{\text{kin}}$ is the initial classical kinetic energy of the trajectory. If the transition is energetically uphill ($\Delta E>E_{\text{kin}}$), the final kinetic energy is negative, and the hop is rejected. The trajectory remains on the active state $c$, which constitutes a frustrated hop.
