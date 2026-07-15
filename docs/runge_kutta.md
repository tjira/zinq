# Runge–Kutta Methods

Runge–Kutta methods are a family of implicit and explicit iterative algorithms used to approximate the solutions of ordinary differential equations. They achieve high-order accuracy by evaluating the derivative of the state vector at multiple intermediate stages within a single time step.

---

## I. Mathematical Formulation

### 1. General Runge–Kutta Steps

For an ordinary differential equation defined by state vector $\mathbf{y}$ and derivative function $\mathbf{f}$ as

$$
\frac{d\mathbf{y}}{dt}=\mathbf{f}(t,\mathbf{y})
$$

a general $s$-stage Runge–Kutta method advances the state from time $t_n$ to $t_{n+1}=t_n+h$ using step size $h$ as

$$
\mathbf{y}_{n+1}=\mathbf{y}_n+h\sum_{i=1}^sb_i\mathbf{k}_i
$$

where $\mathbf{k}_i$ represents the stage derivative vector for stage $i$, computed recursively as

$$
\mathbf{k}_i=\mathbf{f}\left(t_n+c_ih,\mathbf{y}_n+h\sum_{j=1}^{s}a_{ij}\mathbf{k}_j\right)
$$

with coefficients $a_{ij}$, $b_i$, and $c_i$ defining the specific integration method. For explicit methods, the stage summation runs only up to $j=1,\dots,i-1$, meaning the current stage derivative depends only on previously computed derivatives.

### 2. Butcher Tableaus

The algebraic coefficients of a Runge–Kutta method are compactly organized in a structured table called a Butcher tableau. The tableau represents the coefficients as

$$
\begin{array}{c|cccc}c_1&a_{11}&a_{12}&\dots&a_{1s}\\c_2&a_{21}&a_{22}&\dots&a_{2s}\\\vdots&\vdots&\vdots&\ddots&\vdots\\c_s&a_{s1}&a_{s2}&\dots&a_{ss}\\\hline&b_1&b_2&\dots&b_s\end{array}
$$

where the matrix element $a_{ij}$ represents the weight of stage $j$ in calculating stage $i$, $b_i$ represents the weight of stage $i$ in calculating the final step update, and $c_i$ represents the fractional time step for stage $i$.

---

## II. Implementation and Methods

### 3. Compile-Time Optimization and Allocation Tricks

In the codebase, Runge–Kutta steppers are implemented as generic types parameterized by a compile-time `ButcherTableau` structure. This allows the compiler to unroll the stage summation loops and optimize away terms where $a_{ij} = 0$ or $b_i = 0$. To minimize memory fragmentation and allocation overhead during dynamics simulations, the stage derivative vectors $\mathbf{k}_i$ are allocated in a single, contiguous block of memory and accessed as slices. During a step, the temporary state is accumulated in a pre-allocated vector `tmp` before invoking the user-provided derivative evaluation callback.

### 4. First-Order Euler Method

The first-order Euler method (RK1) is the simplest Runge–Kutta scheme, containing a single stage. Its Butcher tableau is given by

$$
\begin{array}{c|c}0&0\\\hline&1\end{array}
$$

which yields a first-order accurate update step.

### 5. Classical Fourth-Order Runge–Kutta Method

The classical fourth-order Runge–Kutta method (RK4) is a widely used four-stage integration scheme. Its Butcher tableau is constructed as

$$
\begin{array}{c|cccc}0&0&0&0&0\\1/2&1/2&0&0&0\\1/2&0&1/2&0&0\\1&0&0&1&0\\\hline&1/6&1/3&1/3&1/6\end{array}
$$

which yields fourth-order global accuracy, providing a robust balance between computational cost and numerical precision.
