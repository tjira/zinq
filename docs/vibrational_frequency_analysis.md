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
