# Dual Numbers and Automatic Differentiation

This document provides a mathematically rigorous explanation of dual numbers and how they enable forward-mode automatic differentiation as implemented in the codebase.

---

## I. Mathematical Definition of Dual Numbers

A dual number $u$ is defined as

$$
u=x+y\epsilon
$$

where $x$ and $y$ are real numbers representing the real and dual parts, respectively, and $\epsilon$ is the dual unit satisfying

$$
\epsilon\neq0,\quad\epsilon^2=0.
$$

In the codebase, dual numbers are implemented via the `ScalarDual(T)` generic structure, where the real part $x$ corresponds to the field `val` and the dual part $y$ corresponds to the field `der`.

### Algebraic Properties and Code Implementation

Given two dual numbers $u=x_1+y_1\epsilon$ and $v=x_2+y_2\epsilon$, the fundamental algebraic operations are mapped directly to methods on the `ScalarDual(T)` structure:

* **Addition**: Evaluated via the `add` method as

$$
u+v=(x_1+x_2)+(y_1+y_2)\epsilon
$$

and scalar addition `adds` as $u+c=(x_1+c)+y_1\epsilon$.

* **Subtraction**: Evaluated via the `sub` method as

$$
u-v=(x_1-x_2)+(y_1-y_2)\epsilon
$$

and scalar subtraction `subs` as $u-c=(x_1-c)+y_1\epsilon$.

* **Multiplication**: Evaluated via the `mul` method using the product rule as

$$
u\cdot v=x_1x_2+(x_1y_2+y_1x_2)\epsilon
$$

and scalar multiplication `muls` as $u\cdot c=x_1c+y_1c\epsilon$.

* **Division**: Evaluated via the `div` method using the quotient rule as

$$
\frac{u}{v}=\frac{x_1}{x_2}+\left(\frac{y_1x_2-x_1y_2}{x_2^2}\right)\epsilon
$$

and scalar division `divs` as $u/c=x_1/c+(y_1/c)\epsilon$, under the condition that the denominator is non-zero.

* **Exponential**: Evaluated via the `exp` method using the chain rule as

$$
e^u=e^{x_1}+e^{x_1}y_1\epsilon
$$

which propagates the derivative through exponential functions.

* **Absolute Value**: Evaluated via the `abs` method as

$$
|u|=|x_1|+\text{sgn}(x_1)y_1\epsilon
$$

which implements the derivative of the absolute value function.

---

## II. Differentiation with Dual Numbers

The connection between dual numbers and differentiation arises from the Taylor series expansion of a real-analytic function $f$ evaluated at a dual number $u=x+y\epsilon$, which is given by

$$
f(x+y\epsilon)=\sum_{k=0}^{\infty}\frac{f^{(k)}(x)}{k!}(y\epsilon)^k=f(x)+f'(x)y\epsilon+\frac{f''(x)}{2!}y^2\epsilon^2+\dots
$$

Since all powers of $\epsilon$ greater than or equal to two are zero, the infinite series truncates exactly after the first-order term, leaving

$$
f(x+y\epsilon)=f(x)+f'(x)y\epsilon.
$$

By setting the dual component $y=1$, evaluating $f(x+\epsilon)$ yields the function value in the real part and the exact derivative in the dual part.

### Generalizing to Multivariate Functions

For a multivariate function $f$ mapping $\mathbb{R}^n$ to $\mathbb{R}$, we can compute the partial derivative with respect to the coordinate $r_j$ by setting the input vector $\mathbf{r}$ to have a dual component of one at index $j$ and zero elsewhere, which we write as

$$
\mathbf{r}_{\text{dual}}=\mathbf{r}+\mathbf{e}_j\epsilon
$$

where $\mathbf{e}_j$ is the $j$-th standard basis vector. Evaluating the function on this dual vector yields

$$
f(\mathbf{r}+\mathbf{e}_j\epsilon)=f(\mathbf{r})+\frac{\partial f(\mathbf{r})}{\partial r_j}\epsilon.
$$

This represents the principle of forward-mode automatic differentiation. The codebase utilizes this trick to evaluate nuclear gradients (such as in Møller–Plesset perturbation theory) by propagating dual numbers throughout the entire self-consistent field and molecular orbital transformation routines, avoiding finite-difference truncation errors.
