# Dual Numbers and Automatic Differentiation

This document provides a highly detailed, mathematically rigorous, yet simple and intuitive explanation of dual numbers and how they enable forward-mode automatic differentiation as implemented in our codebase. If you have ever struggled to understand how a computer can calculate exact derivatives of complex functions without doing symbolic algebra or running into numerical precision errors, this guide is written step-by-step for you.

---

## I. Mathematical Definition of Dual Numbers

To understand dual numbers, it is helpful to think of them as being similar to complex numbers. Just as a complex number has a real part and an imaginary part, a dual number $u$ consists of a real part and a dual part, which we write as

$$
u=x+y\epsilon
$$

where $x$ represents the ordinary real number part, $y$ represents the dual coefficient part, and the symbol $\epsilon$ is a special mathematical unit called the dual unit. The defining rule of the dual unit is that when we square it, it becomes exactly zero, which we write as

$$
\epsilon^2=0
$$

even though the dual unit itself is not equal to zero. In our Zig implementation, we represent these dual numbers using a generic structure called `ScalarDual(T)`. In this structure, the ordinary real part $x$ is stored in the field named `val`, which is short for value, and the dual part $y$ is stored in the field named `der`, which is short for derivative.

### Algebraic Properties and Code Implementation

Because the square of the dual unit is zero, the algebra of dual numbers is simple. If we have two dual numbers, which we denote as $u=x+y\epsilon$ and $v=w+z\epsilon$, we can perform standard arithmetic operations by grouping terms and applying our defining rule to simplify the results. These basic algebraic operations are implemented as methods on the `ScalarDual(T)` structure.

Addition of two dual numbers is performed by adding their real parts and their dual parts separately. We write this sum as

$$
u+v=(x+w)+(y+z)\epsilon
$$

which is implemented in the `add` method. When we add a real scalar constant $c$ to a dual number, we only add it to the real part, which we write as $u+c=(x+c)+y\epsilon$ and implement in the `adds` method.

Subtraction of two dual numbers is performed by subtracting their real parts and their dual parts separately. We write this difference as

$$
u-v=(x-w)+(y-z)\epsilon
$$

which is implemented in the `sub` method. When we subtract a real scalar constant $c$ from a dual number, we only subtract it from the real part, which we write as $u-c=(x-c)+y\epsilon$ and implement in the `subs` method.

Multiplication of two dual numbers is performed by distributing the terms just like we do in ordinary algebra. We expand the product as

$$
u\cdot v=(x+y\epsilon)(w+z\epsilon)=xw+(xz+yw)\epsilon+yz\epsilon^2
$$

where the last term disappears because the square of the dual unit is zero. This simplifies the product to

$$
u\cdot v=xw+(xz+yw)\epsilon
$$

which is implemented in the `mul` method. Notice that the dual part of the product is exactly the product rule from calculus. When we multiply a dual number by a real scalar constant $c$, we multiply both parts, which we write as $u\cdot c=xc+yc\epsilon$ and implement in the `muls` method.

Division of two dual numbers is performed by using a trick similar to rationalizing the denominator in complex numbers. We multiply the numerator and denominator by the conjugate of the denominator, which is $w-z\epsilon$, and simplify the expression as

$$
\frac{u}{v}=\frac{x+y\epsilon}{w+z\epsilon}=\frac{(x+y\epsilon)(w-z\epsilon)}{(w+z\epsilon)(w-z\epsilon)}=\frac{xw+(yw-xz)\epsilon-yz\epsilon^2}{w^2-z^2\epsilon^2}
$$

where the terms containing the squared dual unit disappear. This leaves the simplified quotient as

$$
\frac{u}{v}=\frac{x}{w}+\left(\frac{yw-xz}{w^2}\right)\epsilon
$$

which is implemented in the `div` method. Notice that the dual part is exactly the quotient rule from calculus. For scalar division where we divide a dual number by a real constant $c$, we divide both parts, which we write as $u/c=x/c+(y/c)\epsilon$ and implement in the `divs` method, provided that the denominator is not zero.

Exponential of a dual number is computed by applying the standard properties of the exponential function. We expand the exponential as

$$
e^u=e^{x+y\epsilon}=e^xe^{y\epsilon}
$$

where the second term can be expanded using a Taylor series as $e^{y\epsilon}=1+y\epsilon+(y\epsilon)^2/2!+\dots$ which simplifies to $1+y\epsilon$ because all higher powers of the dual unit are zero. This gives the final result as

$$
e^u=e^x+e^xy\epsilon
$$

which is implemented in the `exp` method. This shows that the derivative of the exponential function propagates naturally.

Absolute value of a dual number is calculated by evaluating the magnitude of the real part and scaling the dual part by the sign of the real part. We write this as

$$
|u|=|x|+\text{sgn}(x)y\epsilon
$$

where $\text{sgn}(x)$ is the sign of the real part, representing the derivative of the absolute value function. This is implemented in the `abs` method.

---

## II. Differentiation with Dual Numbers

The reason dual numbers are so useful for automatic differentiation is because of how they interact with functions. If we have a smooth function $f$ and we evaluate it at a dual number $u=x+y\epsilon$, we can expand the function using a Taylor series around the real part $x$ as

$$
f(x+y\epsilon)=\sum_{k=0}^{\infty}\frac{f^{(k)}(x)}{k!}(y\epsilon)^k=f(x)+f'(x)y\epsilon+\frac{f''(x)}{2!}y^2\epsilon^2+\dots
$$

where $f^{(k)}(x)$ represents the $k$-th derivative of the function evaluated at $x$. Because the square and all higher powers of the dual unit are zero, the infinite sum truncates after the first two terms. This leaves the expression as

$$
f(x+y\epsilon)=f(x)+f'(x)y\epsilon
$$

which means that evaluating the function on a dual number automatically computes both the function value and its exact derivative simultaneously. If we want to find the derivative of a function $f$ at a point $x$, we simply set the input to be the dual number $x+\epsilon$, which has a real part of $x$ and a dual part of one. After running the function, the real part of the result is the function value $f(x)$ and the dual part is the exact derivative $f'(x)$. This process does not suffer from the numerical subtraction errors of finite-difference approximations because no division by a small step size is performed.

### Generalizing to Multivariate Functions

We can extend this technique to compute the derivatives of a function that depends on multiple variables. For a function $f$ that takes a vector of inputs $\mathbf{r}$ and returns a single value, we can compute the partial derivative of the function with respect to the coordinate component $r_j$ by adding the dual unit to that component only. We write the dual input vector as

$$
\mathbf{r}_{\text{dual}}=\mathbf{r}+\mathbf{e}_j\epsilon
$$

where $\mathbf{e}_j$ is a unit vector that has a one at the index $j$ and zeros everywhere else. Evaluating the function on this dual vector yields

$$
f(\mathbf{r}+\mathbf{e}_j\epsilon)=f(\mathbf{r})+\frac{\partial f(\mathbf{r})}{\partial r_j}\epsilon
$$

where the dual part of the output contains the exact partial derivative. This method is called forward-mode automatic differentiation. In our codebase, we use this dual number approach to calculate nuclear gradients for advanced electronic structure methods, such as Møller–Plesset perturbation theory. By running the self-consistent field and molecular orbital transformation routines using the `ScalarDual(T)` type instead of standard real numbers, the program automatically propagates the derivatives throughout the entire calculation, yielding exact analytical gradients without requiring complex manual derivations of the derivative equations.
