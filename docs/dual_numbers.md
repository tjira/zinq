# Dual Numbers and Automatic Differentiation

This document provides a mathematically rigorous explanation of dual numbers, how they enable forward-mode automatic differentiation, and how they are used within the `zinq` library to calculate adiabatic potential energy gradients and Møller–Plesset or Configuration Interaction energy gradients.

---

## 1. Mathematical Definition of Dual Numbers

Dual numbers are an extension of the real numbers, analogous to complex numbers, introduced by William Clifford in 1873.

A dual number $u$ is defined as
$$
u = x + y\epsilon
$$
where $x$ and $y$ are real numbers representing the real and dual parts, respectively, and $\epsilon$ is the dual unit satisfying
$$
\epsilon \neq 0, \quad \epsilon^2 = 0.
$$

### Algebraic Properties

Given two dual numbers $u = x_1 + y_1\epsilon$ and $v = x_2 + y_2\epsilon$, the fundamental algebraic operations are defined by direct expansion.

For addition and subtraction, we have
$$
u \pm v = (x_1 \pm x_2) + (y_1 \pm y_2)\epsilon.
$$

For multiplication, expanding the product yields
$$
u \cdot v = (x_1 + y_1\epsilon)(x_2 + y_2\epsilon) = x_1 x_2 + (x_1 y_2 + y_1 x_2)\epsilon + y_1 y_2\epsilon^2.
$$
Since the dual unit squares to zero, the term containing $\epsilon^2$ vanishes, which simplifies the product to
$$
u \cdot v = x_1 x_2 + (x_1 y_2 + y_1 x_2)\epsilon.
$$

For division, we multiply the numerator and denominator by the conjugate $x_2 - y_2\epsilon$ to obtain
$$
\frac{u}{v} = \frac{(x_1 + y_1\epsilon)(x_2 - y_2\epsilon)}{(x_2 + y_2\epsilon)(x_2 - y_2\epsilon)} = \frac{x_1 x_2 + (y_1 x_2 - x_1 y_2)\epsilon}{x_2^2},
$$
which yields
$$
\frac{u}{v} = \frac{x_1}{x_2} + \left(\frac{y_1 x_2 - x_1 y_2}{x_2^2}\right)\epsilon
$$
under the condition that $x_2 \neq 0$.

---

## 2. Differentiation with Dual Numbers

The connection between dual numbers and differentiation arises from the Taylor series expansion of a real-analytic function $f$ evaluated at a dual number $u = x + y\epsilon$, which is given by
$$
f(x + y\epsilon) = \sum_{k=0}^{\infty} \frac{f^{(k)}(x)}{k!} (y\epsilon)^k = f(x) + f'(x)y\epsilon + \frac{f''(x)}{2!} y^2\epsilon^2 + \dots
$$
Since all powers of $\epsilon$ greater than or equal to two are zero, the infinite series truncates exactly after the first-order term, leaving
$$
f(x + y\epsilon) = f(x) + f'(x)y\epsilon.
$$
By setting the dual component $y = 1$, we obtain
$$
f(x + \epsilon) = f(x) + f'(x)\epsilon.
$$
Thus, evaluating an analytical function over a dual number computes the function value in the real part and the exact derivative in the dual part.

### Generalizing to Multivariate Functions

For a multivariate function $f$ mapping $\mathbb{R}^n$ to $\mathbb{R}$, we can compute the partial derivative with respect to the coordinate $r_j$ by setting the input vector $\mathbf{r}$ to have a dual component of one at index $j$ and zero elsewhere, which we write as
$$
\mathbf{r}_{\text{dual}} = \mathbf{r} + \mathbf{e}_j\epsilon
$$
where $\mathbf{e}_j$ is the $j$-th standard basis vector. Evaluating the function on this dual vector yields
$$
f(\mathbf{r} + \mathbf{e}_j\epsilon) = f(\mathbf{r}) + \frac{\partial f(\mathbf{r})}{\partial r_j}\epsilon.
$$
This represents the principle of forward-mode automatic differentiation. It avoids the truncation errors of numerical finite differences and the expression bloat of symbolic differentiation.
