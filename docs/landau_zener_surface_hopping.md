# Landau–Zener Surface Hopping

Landau–Zener Surface Hopping is a trajectory-based mixed quantum-classical molecular dynamics method that describes nonadiabatic transitions between electronic states at avoided crossings. Unlike Tully's Fewest Switches Surface Hopping which requires explicit nonadiabatic coupling vectors, the Landau–Zener model calculates transition probabilities locally using only the potential energy difference and its curvature along the classical nuclear trajectory.

---

## I. Mathematical Formulation

### 1. Classical Landau–Zener Model

Near an avoided crossing between two adiabatic states, the electronic potential energy difference $Z(t)=|E_j(t)-E_i(t)|$ passes through a local minimum. The classical Landau–Zener model is formulated using a two-state diabatic Hamiltonian where the diabatic energies vary linearly in time as $H_{11}(t)-H_{22}(t)=\alpha t$ and the diabatic coupling $H_{12}$ is constant. The resulting adiabatic energy difference is

$$
Z(t)=\sqrt{\alpha^2t^2+4H_{12}^2}
$$

where the minimum energy gap $Z_{\text{min}}$ occurs at $t=0$, which is

$$
Z_{\text{min}}=2H_{12}
$$

and the second derivative of the energy difference at the crossing point is

$$
\ddot{Z}(0)=\frac{\alpha^2}{2H_{12}}=\frac{\alpha^2}{Z_{\text{min}}}
$$

which yields the linear slope as

$$
\alpha=\sqrt{Z_{\text{min}}\ddot{Z}(0)}
$$

The probability of undergoing a nonadiabatic transition is given by the Landau–Zener transition formula

$$
P=\exp\left(-\frac{2\pi H_{12}^2}{\alpha}\right)
$$

which, by substituting the relations for $H_{12}$ and $\alpha$, can be written entirely in terms of the adiabatic energy difference minimum and its curvature as

$$
P=\exp\left(-\frac{\pi}{2}\sqrt{\frac{Z_{\text{min}}^3}{\ddot{Z}}}\right)
$$

where $Z_{\text{min}}$ is the minimum energy gap at the avoided crossing and $\ddot{Z}$ is the second derivative of the energy difference at that minimum.

---

## II. Implementation and Crossing Detection

### 2. History Buffer and Avoided Crossings

To evaluate transitions, the codebase implements a history buffer trick where the electronic Hamiltonian (or eigenvalues) at the last three time steps is kept in memory. This allows the algorithm to detect avoided crossings and crossings without needing analytical derivatives.

### 3. Adiabatic Transition Probability

In the adiabatic representation, the energy gap $Z(t) = |E_j(t) - E_c(t)|$ is monitored over the last three steps $Z_0$, $Z_1$, and $Z_2$. An avoided crossing is detected when $Z_1 < Z_0$ and $Z_2 \ge Z_1$, indicating a local minimum. The first and second derivatives of the gap are approximated using finite differences, and a parabola is fitted to find the minimum gap $Z_{\text{min}}$ and its curvature $\ddot{Z}$ as

$$
Z_{\text{min}}=Z_1-\frac{b^2}{4a}
$$

where $a = \ddot{Z}/2$ and $b = (Z_2 - Z_0) / (2dt)$, which are then substituted into the exponential Landau–Zener formula to compute the transition probability.

### 4. Diabatic Transition Probability

In the diabatic representation, the energy difference $G(t) = H_{cc}(t) - H_{jj}(t)$ between the active state $c$ and target state $j$ is monitored. A crossing is detected when the sign of $G(t)$ changes between consecutive steps ($G_1 \cdot G_2 \le 0$). The rate of change of the gap is evaluated as $d_G = |G_2 - G_1| / dt$, and the transition probability is calculated using the diabatic coupling element $V_{cj} = |H_{cj}|$ at the crossing point as

$$
P=1-\exp\left(-\frac{2\pi V_{cj}^2}{d_G}\right)
$$

which governs diabatic state transitions.

### 5. Momentum Rescaling and Energy Conservation

To conserve total energy during a nonadiabatic transition, the classical nuclear kinetic energy must compensate for the change in electronic potential energy. When a trajectory hops from state $c$ to state $j$, the change in potential energy is $\Delta E=E_j-E_c$. The nuclear momentum $\mathbf{p}$ is rescaled isotropically by a scaling factor

$$
\gamma=\sqrt{\frac{E_{\text{kin}}-\Delta E}{E_{\text{kin}}}}
$$

where $E_{\text{kin}}$ is the initial classical kinetic energy of the trajectory. If the transition is energetically uphill ($\Delta E>E_{\text{kin}}$), the final kinetic energy is negative, and the hop is rejected. The trajectory remains on the active state $c$, which constitutes a frustrated hop.
