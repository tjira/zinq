# Møller–Plesset Perturbation Theory

Møller–Plesset perturbation theory (MPPT) is a post-Hartree–Fock method that improves upon the mean-field approximation by adding electron correlation effects through Rayleigh–Schrödinger perturbation theory. In this approach, the unperturbed Hamiltonian $\hat{H}_0$ is defined as the sum of the one-electron Fock operators, and the electron correlation is treated as a perturbation.

---

## I. Perturbation Expansion and Spin Variants

### 1. Rayleigh–Schrödinger Perturbation Theory

The total Hamiltonian $\hat{H}$ is partitioned into an unperturbed part $\hat{H}_0$ and a perturbation $\hat{V}$ as

$$
\hat{H}=\hat{H}_0+\hat{V}
$$

where the reference state $| \Phi_0 \rangle$ is an eigenstate of $\hat{H}_0$ with eigenvalue $E_0^{(0)}$ and the unperturbed energy of any determinant $i$ is the sum of its occupied spin-orbital energies $\epsilon_p$ as

$$
E_i^{(0)}=\sum_{p\in\text{det}_i}\epsilon_p
$$

which defines the zeroth-order energy spectrum of the system in terms of the occupied spin-orbital energies $\epsilon_p$ of orbital $p$. The Rayleigh–Schrödinger perturbation theory expands the exact electronic energy and wavefunction in powers of a coupling parameter, leading to iterative corrections at each order. Under the intermediate normalization scheme where $\langle \Phi_0 | \Psi^{(k)} \rangle = \delta_{k0}$, the energy correction of order $k$ is given by

$$
E^{(k)}=\sum_{j\neq0}V_{0j}C_j^{(k-1)}
$$

for $k\ge2$, with the first-order correction being $E^{(1)}=V_{00}$, where $V_{ij}=\langle\Phi_i|\hat{V}|\Phi_j\rangle$ are the perturbation matrix elements in the Slater determinant basis representing the fluctuation potential. The corresponding first-order wavefunction coefficients for $i\neq0$ are given by

$$
C_i^{(1)}=\frac{V_{i0}}{E_0^{(0)}-E_i^{(0)}}
$$

and the coefficients for higher orders $k\ge2$ are updated recursively as

$$
C_i^{(k)}=\frac{\sum_{j\neq0}V_{ij}C_j^{(k-1)}-\sum_{j=1}^{k-1}E^{(j)}C_i^{(k-j)}}{E_0^{(0)}-E_i^{(0)}}
$$

which systematically generates the perturbation expansion, where the energy denominator $E_0^{(0)}-E_i^{(0)}$ represents the excitation energy of determinant $i$ from the reference state. In the codebase, this recurrence is implemented as a general order-$k$ solver that evaluates the coefficients and energy corrections recursively up to any specified perturbation order.

### 2. Generalized and Spin-Unrestricted Møller–Plesset Variants

The formulation of Møller–Plesset perturbation theory depends on the choice of the reference wavefunction. The codebase supports two spin formulations: Restricted Møller–Plesset (RMP) theory, which uses a closed-shell Restricted Hartree–Fock reference, and Generalized Møller–Plesset (GMP) theory, which uses a Generalized Hartree–Fock reference. GMP accommodates non-collinear spin and spin-orbit coupling, and is formulated entirely in a general spin-orbital basis where the spin-blocking structure is relaxed.

---

## II. Nuclear Derivatives

### 3. Analytical Nuclear Gradient

The analytical gradient of the correlation energy of order $k$ with respect to a nuclear coordinate $x$ is formulated by propagating the exact derivatives of the molecular orbital coefficients, orbital energies, and molecular integrals through the perturbation equations. Since the Møller–Plesset wavefunction is not variationally optimized with respect to the molecular orbital coefficients, the nuclear gradient of the correlation energy depends explicitly on the nuclear derivative of the molecular orbital coefficients. By differentiating the energy expression of order $k$, the gradient is given by

$$
\frac{dE^{(k)}}{dx}=\sum_{j\neq0}\left(\frac{dV_{0j}}{dx}C_j^{(k-1)}+V_{0j}\frac{dC_j^{(k-1)}}{dx}\right)
$$

where the derivatives of the wavefunction coefficients $dC_j^{(k-1)}/dx$ are evaluated recursively. The orbital energy derivatives and molecular orbital coefficient derivatives are obtained analytically by solving the coupled-perturbed Hartree–Fock (CPHF) equations, and the derivatives of the molecular integrals are evaluated in the atomic orbital basis.

### 4. Analytical Gradients via Dual Numbers

An alternative and highly efficient method to compute the analytical nuclear gradient of the correlation energy is through the use of dual numbers, which form an algebraic system that automates first-derivative calculations. Dual numbers extend real numbers by introducing an infinitesimal unit $\epsilon$ that satisfies $\epsilon^2=0$, so that any differentiable function $f$ evaluated on a dual number yields its value and derivative as

$$
f(x+dx\epsilon)=f(x)+f'(x)dx\epsilon
$$

without finite-difference errors. To apply this to correlation energy gradients, the nuclear coordinates are perturbed as $x\to x+\epsilon$ and the atomic orbital integrals are computed as dual numbers containing their corresponding derivative values. By propagating these dual integrals through the self-consistent field iterations and the molecular orbital transformation, we obtain dual-valued molecular orbital coefficients $\mathbf{C}+\mathbf{C}'\epsilon$ and orbital energies $\epsilon_p+\epsilon'_p\epsilon$. Evaluating the correlation energy expression of order $k$ using these dual-valued quantities directly yields the analytical gradient as the coefficient of the dual unit $\epsilon$, which bypasses the tedious manual formulation of CPHF terms and the explicit derivatives of the molecular integrals. In the codebase, this is executed by performing the entire calculation using the `ScalarDual(T)` generic type, propagating values and first derivatives simultaneously through all SCF and MO transformation loops.
