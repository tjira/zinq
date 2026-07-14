# Møller–Plesset Perturbation Theory

Møller–Plesset perturbation theory (MPPT) is a post-Hartree–Fock method that improves upon the mean-field approximation by adding electron correlation effects through Rayleigh–Schrödinger perturbation theory. In this approach, the unperturbed Hamiltonian $\hat{H}_0$ is defined as the sum of the one-electron Fock operators, and the electron correlation is treated as a perturbation.

---

## 1. Rayleigh–Schrödinger Perturbation Theory

The total Hamiltonian $\hat{H}$ is divided into the unperturbed Hamiltonian $\hat{H}_0$ and the perturbation $\hat{V}$ as

$$
\hat{H}=\hat{H}_0+\hat{V}
$$

where the unperturbed electronic energy of determinant $i$ is the sum of the occupied orbital energies

$$
E_i^{(0)}=\sum_{p\in\text{det}_i}\epsilon_p
$$

with $\epsilon_p$ representing the spin-orbital energy of orbital $p$. The energy corrections at arbitrary order $k$ are calculated iteratively using the perturbation coefficients as

$$
E^{(k)}=\sum_{j\neq0}H_{0j}C_j^{(k-1)}
$$

where $H_{0j}=\langle\text{det}_0|\hat{H}|\text{det}_j\rangle$ are the Hamiltonian matrix elements in the Slater determinant basis, and the perturbation coefficients are updated as

$$
C_i^{(k)}=\frac{\sum_jH_{ij}C_j^{(k-1)}-E_i^{(0)}C_i^{(k-1)}-\sum_{j=1}^kE^{(j)}C_i^{(k-j)}}{E_0^{(0)}-E_i^{(0)}}
$$

which yields the Rayleigh–Schrödinger perturbation series expansion.

---

## 2. Second-Order Møller–Plesset Energy

For the second-order correction, the summation over determinants simplifies to double excitations from occupied orbitals $i,j$ to virtual orbitals $a,b$ due to Brillouin's theorem. The resulting MP2 correlation energy correction is given by

$$
E^{(2)}=\sum_{i<j}^{\text{occ}}\sum_{a<b}^{\text{vir}}\frac{\left|\langle\phi_i\phi_j||\phi_a\phi_b\rangle\right|^2}{\epsilon_i+\epsilon_j-\epsilon_a-\epsilon_b}
$$

where $\langle\phi_i\phi_j||\phi_a\phi_b\rangle$ are the antisymmetrized two-electron integrals in the molecular spin-orbital basis, and the denominator represents the orbital energy differences.

---

## 3. Analytical Nuclear Gradient

The analytical gradient of the Møller–Plesset correlation energy correction with respect to a nuclear coordinate $x$ is formulated by propagating the exact derivatives of the molecular orbital coefficients, orbital energies, and molecular integrals through the perturbation equations. Since the Møller–Plesset wavefunction is not variational, the nuclear gradient of the correlation energy depends explicitly on the nuclear derivative of the molecular orbital coefficients as

$$
\frac{dE^{(2)}}{dx}=\sum_{i<j}^{\text{occ}}\sum_{a<b}^{\text{vir}}\left(\frac{2\langle\phi_i\phi_j||\phi_a\phi_b\rangle}{\epsilon_i+\epsilon_j-\epsilon_a-\epsilon_b}\frac{d\langle\phi_i\phi_j||\phi_a\phi_b\rangle}{dx}-\frac{\left|\langle\phi_i\phi_j||\phi_a\phi_b\rangle\right|^2}{(\epsilon_i+\epsilon_j-\epsilon_a-\epsilon_b)^2}\left(\frac{d\epsilon_i}{dx}+\frac{d\epsilon_j}{dx}-\frac{d\epsilon_a}{dx}-\frac{d\epsilon_b}{dx}\right)\right)
$$

where the orbital energy derivatives $d\epsilon_i/dx$ and molecular orbital coefficient derivatives $dC_{\mu i}/dx$ are obtained analytically by solving the coupled-perturbed Hartree–Fock (CPHF) equations, and the derivatives of the molecular integrals are evaluated in the atomic orbital basis.
