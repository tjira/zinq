# Classical Molecular Dynamics

Classical molecular dynamics (MD) propagates the classical trajectories of atomic nuclei on potential energy surfaces defined by electronic structure calculations. By solving Newton's equations of motion, MD simulates the time-dependent physical movements of atoms to explore conformational spaces and thermodynamic properties.

---

## I. Trajectory Propagation

### 1. Velocity Verlet Integration

The classical degrees of freedom are propagated in time using the Velocity Verlet algorithm, which integrates the coordinates and momenta with second-order accuracy. For coordinates $\mathbf{r}$ and momenta $\mathbf{p}$ of a system with mass $M$, the coordinates are updated as

$$
\mathbf{r}(t+dt)=\mathbf{r}(t)+\frac{\mathbf{p}(t)}{M}dt+\frac{\mathbf{f}(t)}{2M}dt^2
$$

where $dt$ is the integration time step and $\mathbf{f}(t)$ represents the force vector. The momenta are then updated in two half-steps as

$$
\mathbf{p}\left(t+\frac{dt}{2}\right)=\mathbf{p}(t)+\frac{1}{2}\mathbf{f}(t)dt
$$

and

$$
\mathbf{p}(t+dt)=\mathbf{p}\left(t+\frac{dt}{2}\right)+\frac{1}{2}\mathbf{f}(t+dt)dt
$$

which ensures symplectic and energy-conserving trajectory propagation.

### 2. Force Evaluation on Potential Surfaces

The forces acting on the nuclei are computed as the negative gradient of the potential energy surface. For a single-state adiabatic simulation, the potential energy $V(\mathbf{r})$ is the electronic energy at the current geometry, yielding the classical force

$$
\mathbf{f}_i=-\frac{\partial V(\mathbf{r})}{\partial\mathbf{r}_i}
$$

for nuclear coordinates $\mathbf{r}_i$ of atom $i$.

---

## II. Non-Adiabatic Dynamics and Implementation

### 3. Ensemble Trajectory Management

In the codebase, classical trajectories are managed in an `Ensemble(T)` structure containing matrices for nuclear coordinates `r`, momenta `p`, and accelerations `a`, alongside a vector `s` storing the active electronic state index for each trajectory. This allows simultaneous propagation of multiple independent classical paths to calculate ensemble-averaged observables such as average kinetic and potential energies.

### 4. Ehrenfest Dynamics and Surface Hopping

For multi-state systems, the electronic wavefunction is represented as a linear combination of electronic states with time-dependent coefficients. In Ehrenfest dynamics, the classical nuclei propagate on a mean-field potential surface obtained by averaging over the active electronic states. The effective potential energy is evaluated using the density matrix elements $\rho_{kl}=c_k^*c_l$ as

$$
V_{\text{eff}}(\mathbf{r},t)=\sum_{k,l}\rho_{kl}V_{kl}(\mathbf{r})
$$

where $V_{kl}$ represents the electronic potential energy and coupling matrix elements, resulting in a smooth interpolation of forces across multiple states. In surface hopping, trajectories propagate on a single state $c$ and stochastically transition to another state $j$. If a transition is accepted, the momentum is rescaled isotropically using the factor

$$
\gamma=\sqrt{\frac{E_{\text{kin}}-\Delta E}{E_{\text{kin}}}}
$$

where $E_{\text{kin}}$ is the initial kinetic energy and $\Delta E = E_j - E_c$ is the potential energy change. If $E_{\text{kin}} < \Delta E$ during an uphill transition, the hop is rejected (frustrated hop) and the trajectory continues on the original state.
