# Geometry Optimization

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of geometry optimization methods and how they are implemented in our framework. If you want to understand how a computer can automatically adjust the positions of atoms in a molecule to find its most stable, lowest-energy shape, this guide is written step-by-step for you.

---

## I. Optimization Methods

To understand geometry optimization, we must first think of the potential energy of a molecule as a hilly landscape, which we call the potential energy surface. The coordinates of the atoms define our position on this landscape, and the total energy of the molecule is the height of the hills. When we build a molecule in a computer, our initial guess for the positions of the atoms will usually place us somewhere on the side of a hill. In nature, molecules want to be in their most stable state, which corresponds to the lowest possible energy at the bottom of a valley (a local minimum). Geometry optimization is the process of calculating the forces on the atoms and moving them step-by-step downhill until we reach the bottom of the valley.

### 1. Steepest Descent Method

The steepest descent method is the simplest optimization algorithm. It works by calculating the gradient of the energy with respect to the coordinates of the atoms (which is the negative of the force acting on the atoms) and moving the coordinates directly in the downhill direction. At each iteration $k$, the coordinates are updated as

$$
\mathbf{x}_{k+1}=\mathbf{x}_k-\alpha\mathbf{g}_k
$$

where $\mathbf{x}_k$ is a vector containing the coordinates of all the atoms, $\mathbf{g}_k$ is the gradient of the energy at those coordinates, and $\alpha$ is a small number representing the step size. This method is highly robust when we are far away from the minimum because it always points downhill. However, as we get closer to the bottom of the valley, it tends to oscillate back and forth across the valley floor, resulting in very slow convergence.

### 2. Broyden–Fletcher–Goldfarb–Shanno Method

To avoid the oscillations of the steepest descent method, we use the Broyden–Fletcher–Goldfarb–Shanno (BFGS) method. BFGS is a quasi-Newton method that uses information about both the slope (gradient) and the curvature of the energy surface. The curvature is represented by the Hessian matrix, which contains the second derivatives of the energy. Instead of calculating the expensive Hessian matrix from scratch at every step, BFGS builds an approximation of the inverse Hessian matrix, which we denote as $\mathbf{H}_k$. The coordinates are updated at each step as

$$
\mathbf{x}_{k+1}=\mathbf{x}_k-\alpha\mathbf{H}_k\mathbf{g}_k
$$

where the inverse Hessian matrix $\mathbf{H}_k$ scales and rotates the step direction to point directly towards the minimum, allowing the optimizer to take large, efficient steps.

---

## II. Hessian Update and Convergence

### 3. BFGS Update Formula

At the start of the optimization, we do not know the curvature of the energy surface, so we initialize the approximate inverse Hessian matrix $\mathbf{H}$ as the identity matrix. As we take steps, we measure the change in our positions and the change in the gradients to learn about the curvature of the landscape. We define the position change vector as $\mathbf{s}_k=\mathbf{x}_{k+1}-\mathbf{x}_k$ and the gradient change vector as $\mathbf{y}_k=\mathbf{g}_{k+1}-\mathbf{g}_k$. The approximate inverse Hessian matrix is updated at each step using the formula

$$
\mathbf{H}_{k+1}=\mathbf{H}_k-\rho_k\left(\mathbf{s}_k\mathbf{u}_k^{\text{T}}+\mathbf{u}_k\mathbf{s}_k^{\text{T}}\right)+c_k\mathbf{s}_k\mathbf{s}_k^{\text{T}}
$$

where the intermediate variables are defined as $\mathbf{u}_k=\mathbf{H}_k\mathbf{y}_k$, $\rho_k=1/(\mathbf{y}_k^{\text{T}}\mathbf{s}_k)$, and $c_k=(\mathbf{y}_k^{\text{T}}\mathbf{u}_k/\mathbf{y}_k^{\text{T}}\mathbf{s}_k+1)/(\mathbf{y}_k^{\text{T}}\mathbf{s}_k)$. This update formula ensures that the approximate inverse Hessian remains symmetric and positive-definite, guaranteeing that every step we take is guaranteed to go downhill.

### 4. Implementation Details and Safeguards

In our codebase, a molecule's coordinates are stored as a 2D matrix of shape `N_atoms x 3`. For the matrix-vector multiplications required by the BFGS algorithm, we temporarily flatten this matrix into a 1D vector of length `3 * N_atoms`. To prevent numerical errors, the code implements safeguards: if the dot product of the coordinate change and gradient change $\mathbf{y}_k^{\text{T}}\mathbf{s}_k$ is negative or smaller than $10^{-12}$, we skip the BFGS update for that step. This prevents division by zero and ensures that the inverse Hessian does not become corrupted due to numerical noise or tiny steps near the minimum.

### 5. Convergence Criteria

To determine when we have reached the bottom of the valley, we monitor the forces acting on the atoms. The optimization is considered converged when the root-mean-square (RMS) and maximum values of the gradient components fall below a user-specified threshold. At this point, the forces on the atoms are virtually zero, meaning the molecule has reached a stable equilibrium geometry.
