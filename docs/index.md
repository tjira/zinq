# Zinq

Zinq is a lightweight, high-performance scientific computing library written entirely from scratch in the Zig programming language. It is designed for molecular electronic structure theory, quantum chemistry simulations, and quantum time-dependent propagation dynamics. By prioritizing mathematical clarity, physical rigor, and optimal memory management, the library provides an accessible, transparent, and highly optimized codebase for researchers and developers working in scientific computing.

Whether you are a beginner looking to understand the core physics of molecular simulations or an expert developer looking for a fast, memory-safe library for quantum dynamics, Zinq is structured to be easy to read and modify. The codebase contains no hidden tricks or black-box components: every algorithm, from simple matrix diagonalization to complex perturbation theories, is implemented from first principles.

## Structure of the Documentation

The documentation is organized into three main sections to help you find the information you need:

* **Mathematical Foundations**: Explains the underlying mathematical tools used throughout the library, such as dual numbers for exact automatic differentiation and Runge–Kutta methods for solving ordinary differential equations.
* **Electronic Structure**: Describes the methods used to calculate the static energy, forces, and properties of electrons in molecules, including Hartree–Fock theory, Density Functional Theory, Møller–Plesset perturbation theory, and Configuration–Interaction methods.
* **Time Evolution and Dynamics**: Explains how we simulate the movements of molecules and quantum wavepackets over time, covering classical molecular dynamics, Ehrenfest dynamics, surface hopping, and grid-based split-operator wavepacket propagation.

[View the API Reference](code/index.html)
