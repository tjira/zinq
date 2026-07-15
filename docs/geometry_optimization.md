# Geometry Optimization

Geometry optimization is a computational procedure used to find the spatial arrangement of a set of atoms that corresponds to a local minimum of the potential energy surface. This minimum energy configuration represents a stable molecular structure, such as a reactant, product, or intermediate.

---

## I. Optimization Methods

### 1. Steepest Descent Method

The steepest descent method is a first-order optimization algorithm that updates the nuclear coordinates in the direction of the negative energy gradient, which corresponds to the classical force acting on the nuclei. The coordinates are updated at iteration $k$ as

$$
\mathbf{x}_{k+1}=\mathbf{x}_k-\alpha\mathbf{g}_k
$$

where $\mathbf{x}_k$ represents the nuclear coordinates, $\mathbf{g}_k$ represents the energy gradient, and $\alpha$ is a pre-defined step size. This method is highly robust for structures far from equilibrium but converges slowly near the minimum due to oscillatory behavior.

### 2. Broyden–Fletcher–Goldfarb–Shanno Method

The Broyden–Fletcher–Goldfarb–Shanno (BFGS) method is a quasi-Newton optimization algorithm that utilizes both the gradient and an approximate inverse Hessian matrix to determine the step direction. The nuclear coordinates are updated as

$$
\mathbf{x}_{k+1}=\mathbf{x}_k-\alpha\mathbf{H}_k\mathbf{g}_k
$$

where $\mathbf{H}_k$ represents the approximate inverse Hessian matrix at iteration $k$. The inverse Hessian is updated at each step to incorporate curvature information from the changes in coordinates and gradients.

---

## II. Hessian Update and Convergence

### 3. BFGS Update Formula

The approximate inverse Hessian matrix is updated at each iteration using coordinate differences $\mathbf{s}_k=\mathbf{x}_{k+1}-\mathbf{x}_k$ and gradient differences $\mathbf{y}_k=\mathbf{g}_{k+1}-\mathbf{g}_k$. The update formula is constructed as

$$
\mathbf{H}_{k+1}=\mathbf{H}_k-\rho_k\left(\mathbf{s}_k\mathbf{u}_k^{\text{T}}+\mathbf{u}_k\mathbf{s}_k^{\text{T}}\right)+c_k\mathbf{s}_k\mathbf{s}_k^{\text{T}}
$$

where the intermediate variables are defined by $\mathbf{u}_k=\mathbf{H}_k\mathbf{y}_k$, $\rho_k=1/(\mathbf{y}_k^{\text{T}}\mathbf{s}_k)$, and $c_k=(\mathbf{y}_k^{\text{T}}\mathbf{u}_k/\mathbf{y}_k^{\text{T}}\mathbf{s}_k+1)/(\mathbf{y}_k^{\text{T}}\mathbf{s}_k)$. This ensures that the approximate inverse Hessian remains symmetric and positive-definite, directing the search towards a minimum.

### 4. Implementation Details and Safeguards

In the codebase, nuclear coordinates are stored in matrices of shape `N_atoms x 3` and flattened to vectors of length `3 * N_atoms` during matrix-vector operations. The inverse Hessian matrix is initialized as the identity matrix $\mathbf{I}$. To prevent numerical instabilities and division by zero, the BFGS update is skipped if the projection product $\mathbf{y}_k^{\text{T}}\mathbf{s}_k$ is negative or falls below a threshold of $10^{-12}$, which ensures the positive-definiteness of the inverse Hessian is maintained even under noisy gradients or small step sizes.

### 5. Convergence Criteria

The optimization progress is monitored by evaluating the forces acting on the nuclei. The optimization is converged when the root-mean-square (RMS) of the gradient matrix elements falls below a user-specified threshold, indicating that the system has reached a stationary point on the potential energy surface.
