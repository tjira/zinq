# Landau–Zener Surface Hopping

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of the Landau–Zener surface hopping method and how it is implemented in our scientific computing framework. If you want to understand how we can calculate quantum transition probabilities at avoided crossings using only energy levels and their curvature along the classical path, this guide is written step-by-step for you.

---

## I. Mathematical Formulation

To understand Landau–Zener surface hopping, we must look at how molecules transition between electronic states. In Tully's Fewest Switches Surface Hopping (FSSH), the program must calculate the derivative coupling vectors or wavepacket overlaps at every single step, which is computationally expensive. The Landau–Zener model is a simplified alternative that only calculates transition probabilities when the molecule passes through an avoided crossing, which is a point along the trajectory where the energy gap between the active state and another state reaches a local minimum. Instead of calculating couplings at every step, the Landau–Zener method fits a mathematical curve to the energy gap as the trajectory passes through the crossing point and uses this curvature to calculate a single hopping probability.

### 1. Classical Landau–Zener Model

The Landau–Zener model is based on a simple two-state system where the diabatic energies cross each other linearly in time as $H_{11}(t)-H_{22}(t)=\alpha t$ and the coupling between the states is a constant value $H_{12}$. In this model, the adiabatic energy gap $Z(t)$ between the two states is calculated as

$$
Z(t)=\sqrt{\alpha^2t^2+4H_{12}^2}
$$

where the minimum energy gap $Z_{\text{min}}$ occurs at time $t=0$ and is equal to

$$
Z_{\text{min}}=2H_{12}
$$

which is twice the diabatic coupling. The second derivative of this energy gap with respect to time at the crossing point is calculated as

$$
\ddot{Z}(0)=\frac{\alpha^2}{2H_{12}}=\frac{\alpha^2}{Z_{\text{min}}}
$$

which we can rearrange to find the rate of change of the diabatic energy difference as

$$
\alpha=\sqrt{Z_{\text{min}}\ddot{Z}(0)}
$$

where $\alpha$ represents the slope. The probability that the system will undergo a transition between the adiabatic states is given by the Landau–Zener formula

$$
P=\exp\left(-\frac{2\pi H_{12}^2}{\alpha}\right)
$$

where, by substituting the formulas for $H_{12}$ and $\alpha$, we can express the probability entirely in terms of the minimum energy gap and its second derivative as

$$
P=\exp\left(-\frac{\pi}{2}\sqrt{\frac{Z_{\text{min}}^3}{\ddot{Z}}}\right)
$$

where $Z_{\text{min}}$ is the minimum energy gap and $\ddot{Z}$ is the curvature of the gap at the avoided crossing.

---

## II. Implementation and Crossing Detection

### 2. History Buffer and Avoided Crossings

To detect when a trajectory is passing through an avoided crossing, our codebase implements a history buffer: the program stores the electronic energies from the last three time steps in memory. This allows the algorithm to detect when the energy gap reaches a minimum or when the diabatic energies cross without needing analytical derivatives.

### 3. Adiabatic Transition Probability

In the adiabatic representation, the program monitors the energy gap $Z(t) = |E_j(t) - E_c(t)|$ between the active state $c$ and target state $j$ over the last three steps, which we denote as $Z_0$, $Z_1$, and $Z_2$. An avoided crossing is detected when $Z_1 < Z_0$ and $Z_2 \ge Z_1$, indicating that the gap reached a minimum at the middle step. The program then fits a parabola to these three points to find the exact minimum gap $Z_{\text{min}}$ and its curvature $\ddot{Z}$. The minimum gap is calculated as

$$
Z_{\text{min}}=Z_1-\frac{b^2}{4a}
$$

where the coefficients are $a = (Z_2 - 2Z_1 + Z_0) / (2dt^2)$ and $b = (Z_2 - Z_0) / (2dt)$, representing the curvature and slope. The minimum gap and curvature are then plugged into the exponential Landau–Zener formula to compute the hopping probability.

### 4. Diabatic Transition Probability

In the diabatic representation, the program monitors the energy difference $G(t) = H_{cc}(t) - H_{jj}(t)$ between the active state $c$ and target state $j$. A crossing is detected when the sign of $G(t)$ changes between consecutive steps, meaning $G_1 \cdot G_2 \le 0$. The rate of change of the gap is calculated as $d_G = |G_2 - G_1| / dt$, and the transition probability is calculated using the diabatic coupling $V_{cj} = |H_{cj}|$ at the crossing point as

$$
P=1-\exp\left(-\frac{2\pi V_{cj}^2}{d_G}\right)
$$

which governs the probability of transitioning to the other diabatic state.

### 5. Momentum Rescaling and Energy Conservation

Just as in other surface hopping methods, when a trajectory hops from state $c$ to state $j$, the nuclear momentum vector $\mathbf{p}$ must be scaled to conserve energy. The scaling factor is calculated as

$$
\gamma=\sqrt{\frac{E_{\text{kin}}-\Delta E}{E_{\text{kin}}}}
$$

where $E_{\text{kin}}$ is the nuclear kinetic energy and $\Delta E = E_j - E_c$ is the potential energy change. If the hop is uphill and the potential energy change is larger than the available kinetic energy, the hop is rejected (frustrated hop) and the trajectory continues on its original state.
