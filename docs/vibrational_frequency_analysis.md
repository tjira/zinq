# Vibrational Frequency Analysis

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of vibrational frequency analysis and how it is implemented in our framework. If you want to understand how a computer can simulate how a molecule vibrates like a set of balls connected by springs, and how we calculate the frequencies of these vibrations, this guide is written step-by-step for you.

---

## I. Mass-Weighting and Projection

To understand vibrational frequency analysis, we can think of a molecule as a collection of atoms connected by chemical bonds. Around its stable equilibrium geometry, we can model these bonds as harmonic springs. When the atoms are displaced from their equilibrium positions, they will vibrate back and forth. A vibrational frequency analysis calculates the specific frequencies at which these vibrations occur (the vibrational spectrum) and the directions in which the atoms move during each vibration (the normal modes). This analysis also tells us if we have found a true stable minimum (where all frequencies are real and positive) or a transition state (which has one imaginary frequency, indicating the top of an energy barrier).

### 1. Mass-Weighted Hessian Matrix

The curvature of the potential energy surface is described by the Hessian matrix, which contains the second derivatives of the potential energy with respect to the coordinates of the atoms. These second derivatives represent the force constants (the stiffness) of the springs between the atoms. However, heavy atoms move slower than light atoms when subjected to the same force. To account for the masses of the different atoms, we convert the Cartesian Hessian matrix $\mathbf{H}$ into the mass-weighted Hessian matrix $\mathbf{H}^{\text{MW}}$. The elements of this matrix are defined as

$$
H^{\text{MW}}_{i\alpha,j\beta}=\frac{H_{i\alpha,j\beta}}{\sqrt{M_iM_j}}
$$

where $i$ and $j$ represent atomic indices, $\alpha$ and $\beta$ represent Cartesian coordinate directions (x, y, or z), and $M_i$ and $M_j$ represent the masses of the corresponding atoms.

### 2. Projection of Translations and Rotations

A free molecule containing $N$ atoms has $3N$ total degrees of freedom. However, not all of these represent vibrations. The molecule as a whole can translate in three independent directions (along the x, y, and z axes) and rotate as a rigid body in three independent directions (or two for linear molecules). These translations and rotations do not stretch or bend any bonds, so they require zero energy and should have vibrational frequencies of exactly zero. In numerical calculations on a computer, small rounding errors can cause these translation and rotation modes to mix with the true vibrations, yielding small, unphysical non-zero frequencies. To prevent this, we construct an orthonormal projection matrix $\mathbf{P}$ and apply it to the mass-weighted Hessian matrix to yield the projected mass-weighted Hessian matrix

$$
\mathbf{H}^{\text{proj}}=\mathbf{P}\mathbf{H}^{\text{MW}}\mathbf{P}
$$

which forces the translational and rotational frequencies to be exactly zero.

---

## II. Normal Modes and Frequencies

### 3. Implementation Tricks and Projector Construction

To construct the projection matrix $\mathbf{P}$ in our codebase, we perform a step-by-step geometric calculation. First, we find the center of mass of the molecule and translate all the coordinates so that the center of mass lies exactly at the origin. Second, we construct the $3\times3$ inertia tensor of the molecule and diagonalize it using a symmetric eigenvalue solver to find the principal axes of rotation and the principal moments of inertia. Third, we construct the three translational and three rotational basis vectors in the $3N$-dimensional mass-weighted coordinate space. We then orthonormalize these vectors using a Gram–Schmidt procedure to build a transformation matrix $\mathbf{U}_{\text{tr}}$ containing the translation-rotation basis. Finally, the projection matrix is evaluated as

$$
P_{ij}=\delta_{ij}-\sum_kU_{\text{tr},ik}U_{\text{tr},jk}
$$

where $\delta_{ij}$ is the Kronecker delta (which is one if $i=j$ and zero otherwise). We apply this projection matrix to the mass-weighted Hessian using matrix multiplication to remove all translational and rotational contamination.

### 4. Normal Mode Diagonalization

We find the vibrational frequencies and normal modes by diagonalizing the projected mass-weighted Hessian matrix, solving the eigenvalue equation

$$
\mathbf{H}^{\text{proj}}\mathbf{Q}_k=\lambda_k\mathbf{Q}_k
$$

where $\mathbf{Q}_k$ is the eigenvector representing the normal mode coordinate vector for mode $k$, and $\lambda_k$ is the corresponding eigenvalue. The harmonic vibrational frequency $\omega_k$ in atomic units is calculated from the eigenvalue as

$$
\omega_k=\text{sgn}(\lambda_k)\sqrt{|\lambda_k|}
$$

where the sign function is used to preserve negative eigenvalues. A negative eigenvalue yields a negative (imaginary) frequency, which physically represents a direction in coordinate space where the energy decreases, indicating that the molecule is at a transition state (the top of a saddle point) rather than a stable minimum. The frequencies are then converted from atomic units to standard spectroscopic units of wavenumbers ($\text{cm}^{-1}$) using physical constants defined in the code, allowing direct comparison with experimental infrared spectra.

---

## III. Thermochemical Analysis

### 5. Ideal Gas and Rigid-Rotor Harmonic-Oscillator (RRHO) Model

Thermochemical analysis connects quantum electronic structure calculations and vibrational frequencies with macroscopic thermodynamic observables such as internal energy, enthalpy, entropy, and Gibbs free energy. Under the standard ideal-gas, rigid-rotor, and harmonic-oscillator (RRHO) approximations, the total molecular partition function factorizes into independent translational, rotational, vibrational, and electronic components as

$$
q_{\text{tot}}=q_{\text{trans}}q_{\text{rot}}q_{\text{vib}}q_{\text{elec}}
$$

which allows each thermodynamic property to be evaluated as a sum of separable contributions at a given temperature $T$ and pressure $P$.

### 6. Translational and Rotational Contributions

The translational partition function per unit volume is derived from the quantum particle-in-a-box model for a molecule of total mass $M$ in the continuum limit as

$$
\frac{q_{\text{trans}}}{V}=\left(\frac{2\pi M k_B T}{h^2}\right)^{3/2}
$$

which yields the classical translational thermal energy $E_{\text{trans}}=\frac{3}{2}k_BT$. Applying the Sackur–Tetrode relation at pressure $P$ gives the translational entropy

$$
S_{\text{tran}}=k_B\left[\ln\left(\frac{q_{\text{trans}}}{V}\frac{k_BT}{P}\right)+\frac{5}{2}\right]
$$

where $k_B$ is the Boltzmann constant and $h$ is Planck's constant.

For molecular rotations, the principal moments of inertia $I_A\le I_B\le I_C$ are obtained by diagonalizing the molecular inertia tensor with respect to the center of mass. For a single atom, rotational contributions vanish identically. For a linear molecule with average moment of inertia $I_{\text{lin}}=\frac{I_B+I_C}{2}$, the rotational thermal energy is $E_{\text{rot}}=k_BT$ and the rotational entropy is

$$
S_{\text{rota}}=k_B(\ln q_{\text{rot}}+1)
$$

with the linear rotational partition function defined as

$$
q_{\text{rot}}=\frac{8\pi^2I_{\text{lin}}k_BT}{h^2}
$$

while for a non-linear polyatomic molecule, the rotational thermal energy is $E_{\text{rot}}=\frac{3}{2}k_BT$ and the rotational entropy is

$$
S_{\text{rota}}=k_B\left(\ln q_{\text{rot}}+\frac{3}{2}\right)
$$

where the non-linear rotational partition function is evaluated as

$$
q_{\text{rot}}=\sqrt{\pi}\left(\frac{8\pi^2k_BT}{h^2}\right)^{3/2}\sqrt{I_AI_BI_C}
$$

incorporating all three principal moments of inertia.

### 7. Vibrational and Electronic Contributions

Every real harmonic vibrational normal mode with wavenumber $\tilde{\nu}_k>1\text{ cm}^{-1}$ has energy $\varepsilon_k=hc\tilde{\nu}_k$. The quantum zero-point vibrational energy (ZPVE) representing the persistent vibrational ground-state energy at absolute zero is

$$
E_{\text{zpve}}=\sum_k\frac{1}{2}\varepsilon_k
$$

and using the dimensionless thermal parameter $x_k=\frac{\varepsilon_k}{k_BT}$, the additional thermal vibrational energy arising from thermal population of excited vibrational levels is

$$
E_{\text{vib}}=\sum_k\frac{\varepsilon_k}{e^{x_k}-1}
$$

with the vibrational entropy given by

$$
S_{\text{vibr}}=k_B\sum_k\left[\frac{x_k}{e^{x_k}-1}-\ln(1-e^{-x_k})\right]
$$

summed over all non-imaginary vibrational modes.

Assuming a non-degenerate ground electronic state with spin multiplicity $M_{\text{spin}}=2S+1$, the electronic entropy accounts for the spin multiplicity degeneracy as

$$
S_{\text{elec}}=k_B\ln M_{\text{spin}}
$$

and summing all components yields the total molecular entropy

$$
S_{\text{tot}}=S_{\text{tran}}+S_{\text{rota}}+S_{\text{vibr}}+S_{\text{elec}}
$$

which is converted to molar entropy units of $\text{cal}/(\text{mol}\cdot\text{K})$ by multiplying by the gas constant $R$ divided by the thermochemical calorie conversion factor.

### 8. Enthalpy, Free Energy, and Total State Functions

The thermal correction to the internal energy $E_{\text{thrm}}$ represents the sum of the zero-point vibrational energy and all thermal motions according to

$$
E_{\text{thrm}}=E_{\text{zpve}}+E_{\text{trans}}+E_{\text{rot}}+E_{\text{vib}}
$$

from which the thermal correction to the enthalpy $H_{\text{corr}}$ is obtained by adding the ideal-gas pressure-volume work term $k_BT$ as

$$
H_{\text{corr}}=E_{\text{thrm}}+k_BT
$$

and the thermal correction to the Gibbs free energy $G_{\text{corr}}$ is obtained by incorporating entropic stabilization as

$$
G_{\text{corr}}=H_{\text{corr}}-TS_{\text{tot}}
$$

Combining these thermal corrections with the electronic energy $E_{\text{elec}}$ computed by an electronic structure method (such as Hartree–Fock, DFT, Møller–Plesset perturbation theory, or Configuration Interaction) provides the final macroscopic thermodynamic state functions

$$
E_0=E_{\text{elec}}+E_{\text{zpve}}
$$

for the zero-point corrected energy,

$$
E_{\text{tot}}=E_{\text{elec}}+E_{\text{thrm}}
$$

for the total internal energy,

$$
H_{\text{tot}}=E_{\text{elec}}+H_{\text{corr}}
$$

for the total enthalpy, and

$$
G_{\text{tot}}=E_{\text{elec}}+G_{\text{corr}}
$$

for the total Gibbs free energy at temperature $T$ and pressure $P$.
