# Vibrational Frequency Analysis

Vibrational frequency analysis is used to determine the normal modes and corresponding vibrational frequencies of a molecular system. This calculation is essential for identifying whether a stationary point is a minimum or a transition state, and for computing thermodynamic properties.

---

## I. Mass-Weighting and Projection

### 1. Mass-Weighted Hessian Matrix

The vibrational frequencies are obtained from the second derivatives of the energy with respect to nuclear coordinates, which form the Hessian matrix. To decouple the motion of different atoms with varying masses, the Cartesian Hessian matrix $\mathbf{H}$ is converted into the mass-weighted Hessian matrix $\mathbf{H}^{\text{MW}}$ with elements defined by

$$
H^{\text{MW}}_{i\alpha,j\beta}=\frac{H_{i\alpha,j\beta}}{\sqrt{M_iM_j}}
$$

where $i$ and $j$ represent atomic indices, $\alpha$ and $\beta$ represent Cartesian coordinate directions, and $M_i$ and $M_j$ represent the corresponding atomic masses.

### 2. Projection of Translations and Rotations

A free molecule has three translational and three rotational (two for linear molecules) degrees of freedom that have zero frequency. To prevent numerical mixing with the true vibrational modes, these degrees of freedom are projected out. An orthonormal projection matrix $\mathbf{P}$ is constructed from center-of-mass and inertia tensor coordinates, and the projected mass-weighted Hessian is calculated as

$$
\mathbf{H}^{\text{proj}}=\mathbf{P}\mathbf{H}^{\text{MW}}\mathbf{P}
$$

which constrains translations and rotations to have zero-valued eigenvalues.

---

## II. Normal Modes and Frequencies

### 3. Implementation Tricks and Projector Construction

To construct the translation-rotation projection matrix $\mathbf{P}$, the codebase first computes the center of mass of the molecule and translates the coordinates. It then constructs the $3\times3$ inertia tensor and diagonalizes it using a symmetric eigenvalue solver to obtain the principal moments of inertia and the corresponding principal axes vectors. The translational and rotational basis vectors are constructed in the mass-weighted Cartesian space and orthonormalized using a Gram–Schmidt procedure to build a transformation matrix $\mathbf{U}_{\text{tr}}$. The projection matrix is then evaluated as

$$
P_{ij}=\delta_{ij}-\sum_kU_{\text{tr},ik}U_{\text{tr},jk}
$$

which is applied to the mass-weighted Hessian via matrix multiplication as $\mathbf{H}^{\text{proj}}=\mathbf{P}\mathbf{H}^{\text{MW}}\mathbf{P}$, ensuring that translational and rotational contamination is projected out.

### 4. Normal Mode Diagonalization

The projected mass-weighted Hessian is diagonalized to obtain its eigenvalues $\lambda_k$ and eigenvectors as

$$
\mathbf{H}^{\text{proj}}\mathbf{Q}_k=\lambda_k\mathbf{Q}_k
$$

where $\mathbf{Q}_k$ represents the normal mode coordinate vector for mode $k$. The harmonic vibrational frequency $\omega_k$ in atomic units is calculated from the eigenvalue as

$$
\omega_k=\text{sgn}(\lambda_k)\sqrt{|\lambda_k|}
$$

where the sign function preserves negative eigenvalues as imaginary frequencies, representing transition states on the potential energy surface. The frequencies are then converted to wavenumbers in $\text{cm}^{-1}$ using physical constants defined in the codebase for direct comparison with experimental spectra.
