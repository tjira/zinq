# Landau–Zener Surface Hopping

Landau–Zener Surface Hopping is a trajectory-based mixed quantum-classical molecular dynamics method that describes nonadiabatic transitions between electronic states at avoided crossings. Unlike Tully's Fewest Switches Surface Hopping which requires explicit nonadiabatic coupling vectors, the Landau–Zener model calculates transition probabilities locally using only the adiabatic potential energy difference and its curvature along the classical nuclear trajectory.

---

## 1. Mathematical Formulation

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

## 2. Detection of Avoided Crossings and Interpolation

During classical trajectory propagation, the energy difference $Z(t)$ between the active state $c$ and any other state $j$ is monitored. An avoided crossing is detected when the first derivative of the energy difference $\dot{Z}(t)$ changes sign from negative to positive, indicating a local minimum. To find the exact minimum energy gap $Z_{\text{min}}$ and its curvature $\ddot{Z}$ from discrete time steps, a quadratic interpolation of the energy difference is performed as

$$
Z_{\text{min}}=Z(t-\Delta t)-\frac{\dot{Z}^2}{2\ddot{Z}}
$$

where the first derivative $\dot{Z}$ and second derivative $\ddot{Z}$ are approximated from the energy differences at consecutive time steps.

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
