# Runge–Kutta Methods

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of Runge–Kutta methods and how they are implemented in our scientific computing framework. If you have ever wondered how a computer can predict the future behavior of a physical system by solving differential equations step-by-step, this guide is designed for you.

---

## I. Mathematical Formulation

To understand Runge–Kutta methods, we must first understand what an ordinary differential equation is. An ordinary differential equation describes how a system changes over time. For example, if we have a state vector $\mathbf{y}$ that represents the positions, velocities, or quantum amplitudes of a system, its rate of change with respect to time $dt$ is given by a derivative function $\mathbf{f}$ as

$$
\frac{d\mathbf{y}}{dt}=\mathbf{f}(t,\mathbf{y})
$$

where $t$ is the time and $\mathbf{y}$ is the state of the system. Because we usually cannot solve these equations exactly by hand, we use a computer to estimate the state of the system at future times. The simplest way to do this is to start at a time $t_n$ with a known state $\mathbf{y}_n$ and take a small step forward in time of size $h$, which we call the time step, to find the new state $\mathbf{y}_{n+1}$ at time $t_{n+1}=t_n+h$.

### 1. General Runge–Kutta Steps

Rather than just calculating the derivative at the beginning of the step and moving in that direction, which is inaccurate, a Runge–Kutta method calculates the derivative at multiple intermediate points, which we call stages. By taking a weighted average of these intermediate derivatives, we get a much more accurate prediction of the final state. A general $s$-stage Runge–Kutta method computes the next state $\mathbf{y}_{n+1}$ from the current state $\mathbf{y}_n$ as

$$
\mathbf{y}_{n+1}=\mathbf{y}_n+h\sum_{i=1}^sb_i\mathbf{k}_i
$$

where $s$ is the total number of stages, $b_i$ represents the weight of stage $i$ in the final update, and $\mathbf{k}_i$ represents the stage derivative vector for stage $i$. Each stage derivative $\mathbf{k}_i$ is computed recursively as

$$
\mathbf{k}_i=\mathbf{f}\left(t_n+c_ih,\mathbf{y}_n+h\sum_{j=1}^{s}a_{ij}\mathbf{k}_j\right)
$$

where $a_{ij}$ represents the weight of stage $j$ in calculating the temporary state for stage $i$, and $c_i$ represents the fraction of the time step $h$ at which the stage is evaluated. In an explicit Runge–Kutta method, the sum for stage $i$ only goes up to $j=1,\dots,i-1$, which means we only use stage derivatives that we have already calculated. This makes the method easy to calculate step-by-step.

### 2. Butcher Tableaus

To make it easy to write computer code that can run any Runge–Kutta method, we organize the coefficients $a_{ij}$, $b_i$, and $c_i$ into a standard table called a Butcher tableau. The tableau is structured as

$$
\begin{array}{c|cccc}c_1&a_{11}&a_{12}&\dots&a_{1s}\\c_2&a_{21}&a_{22}&\dots&a_{2s}\\\vdots&\vdots&\vdots&\ddots&\vdots\\c_s&a_{s1}&a_{s2}&\dots&a_{ss}\\\hline&b_1&b_2&\dots&b_s\end{array}
$$

where the left column contains the time fractions $c_i$, the main matrix contains the stage weights $a_{ij}$, and the bottom row contains the final weights $b_i$. If a method is explicit, all the matrix entries on and above the main diagonal are zero because a stage cannot depend on its own derivative or future stage derivatives.

---

## II. Implementation and Methods

### 3. Compile-Time Optimization and Allocation Tricks

In our codebase, we want the numerical integration to run as fast as possible. To achieve this, we implement the Runge–Kutta steppers as generic types that take a compile-time `ButcherTableau` structure as a parameter. By using Zig's compile-time features, the compiler can inspect the coefficients of the tableau before generating the machine code. If a coefficient $a_{ij}$ or $b_i$ is equal to zero, the compiler completely removes that term from the calculation, avoiding unnecessary multiplications and memory lookups. Furthermore, to avoid allocating and deallocating memory during a simulation, which slows down the program, all stage derivative vectors $\mathbf{k}_i$ are allocated once in a single contiguous block of memory when the simulation starts. The program uses small slices of this pre-allocated memory to store each stage derivative, and a single temporary state vector `tmp` is used to hold intermediate calculations before calling the user's derivative function.

### 4. First-Order Euler Method

The first-order Euler method (RK1) is the simplest possible Runge–Kutta method. It has only one stage, and its Butcher tableau is written as

$$
\begin{array}{c|c}0&0\\\hline&1\end{array}
$$

which yields a simple step update where we evaluate the derivative at the start of the step and assume it remains constant across the entire interval. This method is simple but has a large error that grows linearly with the time step size.

### 5. Classical Fourth-Order Runge–Kutta Method

The classical fourth-order Runge–Kutta method (RK4) is one of the most popular integration schemes in scientific computing because it is highly accurate without being too expensive. It uses four stages, and its Butcher tableau is written as

$$
\begin{array}{c|cccc}0&0&0&0&0\\1/2&1/2&0&0&0\\1/2&0&1/2&0&0\\1&0&0&1&0\\\hline&1/6&1/3&1/3&1/6\end{array}
$$

which translates to the following step-by-step algorithm. First, we compute the derivative $\mathbf{k}_1$ at the beginning of the step. Second, we use $\mathbf{k}_1$ to estimate the state of the system at the middle of the time step and compute the derivative $\mathbf{k}_2$ there. Third, we use $\mathbf{k}_2$ to make a better estimate of the state at the middle of the time step and compute the derivative $\mathbf{k}_3$. Fourth, we use $\mathbf{k}_3$ to estimate the state at the end of the time step and compute the derivative $\mathbf{k}_4$. Finally, we combine these four derivatives as a weighted average to calculate the new state. This method has an error that decreases with the fourth power of the time step size, making it extremely precise for molecular and quantum dynamics.
