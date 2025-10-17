---
title: Hilbert Spaces and Quantum States
parent: Mathematical Background
layout: default
nav_order: 2
---
{% include mathjax.html %}

# Hilbert Spaces and Quantum States<!--\label{sec:hilbert_spaces}-->

Hilbert space is a fundamental concept in quantum mechanics, providing the mathematical framework for describing quantum states and their evolution. In this section, we will explore the definition of Hilbert spaces, their properties, and how they relate to quantum states. We will start with the definition.

{:.definition}
> A Hilbert space $$\mathcal{H}$$ is a complete inner product space, which means it is a vector space equipped with an inner product that allows for the definition of length and angle, and it is complete with respect to the norm induced by the inner product.

Some examples of Hilbert spaces include the finite-dimensional vector spaces $$\mathbb{R}^n$$ or $$\mathbb{C}^n$$ or, as is often the case in quantum mechanics, the space of all square-integrable functions $$L^2(\mathbb{R}^n)$$. Formally

$$
\begin{equation}
L^2(\mathbb{R}^n)=\left\lbrace\psi :\mathbb{R}^n\to\mathbb{C}\middle|\int_{\mathbb{R}^n}\psi^*(\mathbf{r})\psi(\mathbf{r})\mathrm{d}\mathbf{r}<\infty\right\rbrace
\end{equation}
$$

where $$\psi^*(x)$$ is the complex conjugate of $$\psi(x)$$. When dealing with quantum states, we often use the Dirac notation, where quantum states are represented as ket vectors $$\ket{\psi}$$. The inner product between two elements $$\ket{\psi}$$ and $$\ket{\phi}$$ in the Hilbert space $$L^2(\mathbb{R}^n)$$ is defined as

$$
\begin{equation}
\braket{\psi\vert\phi}=\int_{\mathbb{R}^n}\psi^*(\mathbf{r})\phi(\mathbf{r})\,\mathrm{d}\mathbf{r}.
\end{equation}
$$
