# Classical Molecular Dynamics

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of classical molecular dynamics (MD) and how it is implemented in our framework. If you want to understand how we can simulate the physical movement of atoms in a molecule over time by solving Newton's equations of motion, this guide is written step-by-step for you.

---

## I. Trajectory Propagation

To understand molecular dynamics, we must look at how we partition a molecule. A molecule is made of light, fast-moving electrons and heavy, slow-moving atomic nuclei. Because the nuclei are thousands of times heavier than the electrons, we can make a simplifying assumption: we can treat the nuclei as classical point-like particles that move according to Newton's laws of motion, while the electrons are treated quantum mechanically. The electrons adjust instantly to any movement of the nuclei, creating an average potential energy surface on which the nuclei move. Classical molecular dynamics is the process of calculating the forces acting on the nuclei at their current positions and moving them forward in time using small steps to simulate their physical trajectories.

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

## II. Non-Adiabatic Dynamics and Implementation

### 3. Ensemble Trajectory Management

In many simulations, running a single molecular trajectory is not enough to get reliable results because molecular behavior is statistical. To calculate observable properties (like average kinetic energy or chemical reaction rates), we run a large group of independent trajectories, which we call an ensemble. In our codebase, these trajectories are managed using an `Ensemble(T)` structure. This structure contains coordinate matrix `r`, momentum matrix `p`, and acceleration matrix `a` for all trajectories, along with a state vector `s` that stores the active electronic state index for each independent path.

### 4. Ehrenfest Dynamics and Surface Hopping

When a molecule is excited by light, the electrons can transition between different energy states. To simulate these non-adiabatic processes where multiple electronic states are involved, we use two main methods. In Ehrenfest dynamics, the electronic wavefunction is represented as a linear combination of electronic states with complex coefficients. The classical nuclei move on an average potential surface obtained by taking a weighted sum over the active states. The effective potential energy is evaluated using the density matrix elements $\rho_{kl}=c_k^*c_l$ as

$$
V_{\text{eff}}(\mathbf{r},t)=\sum_{k,l}\rho_{kl}V_{kl}(\mathbf{r})
$$

where $V_{kl}$ represents the electronic energy and coupling elements, resulting in a smooth average force. In surface hopping, each trajectory propagates on a single active electronic state $c$ and can undergo sudden, stochastic transitions to another state $j$. If a transition is accepted, the nuclear momentum $\mathbf{p}$ must be rescaled to conserve total energy. The scaling factor is calculated as

$$
\gamma=\sqrt{\frac{E_{\text{kin}}-\Delta E}{E_{\text{kin}}}}
$$

where $E_{\text{kin}}$ is the initial kinetic energy and $\Delta E = E_j - E_c$ is the change in potential energy. If the transition is energetically uphill and the kinetic energy is smaller than the potential energy change, the hop is rejected. The trajectory continues on its original state, which constitutes a frustrated hop.
