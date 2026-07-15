# Zinq

Zinq is a lightweight, high-performance Zig framework designed for molecular electronic structure theory, quantum chemistry simulations, and quantum time propagation dynamics. Built entirely from scratch, the library prioritizes mathematical clarity, physical rigor, and optimal memory management, offering an accessible and transparent codebase for researchers and developers in scientific computing.

[View the API Reference](code/index.html)

---

## Documentation Sections

The **[Mathematical Foundations](dual_numbers.md)** section details the numerical machinery of the library, featuring [Dual Numbers](dual_numbers.md) for automatic differentiation and [Runge–Kutta Methods](runge_kutta.md) for ordinary differential equation integration.

The **[Electronic Structure](density_functional_theory.md)** section covers the quantum chemistry methods, including [Density Functional Theory](density_functional_theory.md), the [Hartree–Fock Method](hartree_fock.md), [Møller–Plesset Perturbation Theory](moller_plesset.md), [Configuration–Interaction](configuration_interaction.md), [Geometry Optimization](geometry_optimization.md), [Vibrational Frequency Analysis](vibrational_frequency_analysis.md), and [Population Analysis](population_analysis.md).

The **[Time Evolution & Dynamics](classical_molecular_dynamics.md)** section describes mixed quantum-classical and grid-based propagation methods, including [Classical Molecular Dynamics](classical_molecular_dynamics.md), [Ehrenfest Dynamics](ehrenfest_dynamics.md), [Fewest Switches Surface Hopping](fewest_switches_surface_hopping.md), [Landau–Zener Surface Hopping](landau_zener_surface_hopping.md), and the [Split-Operator Method](split_operator_method.md).

---

## Key Features

### Mathematical Foundations

The library provides essential mathematical tools tailored for scientific computing. This includes a generic dual number system for exact forward-mode automatic differentiation, bypassing the truncation errors of finite differences. It also features a generic ordinary differential equation solver parameterized by compile-time Butcher tableaus to run optimized Runge–Kutta integration steps.

### Electronic Structure Methods

The electronic structure suite contains Hartree–Fock and Density Functional Theory methods. It supports Restricted and Generalized spin variants, utilizing a molecular grid built from radial shells and Lebedev angular spheres to numerically integrate exchange-correlation functionals. The post-Hartree–Fock correlation models include Configuration–Interaction and arbitrary-order Møller–Plesset perturbation theory. The suite also provides BFGS and steepest descent geometry optimizers, mass-weighted vibrational frequency analysis, and Mulliken charge population analysis.

### Quantum Dynamics and Time Evolution

The dynamics suite supports quantum wavepacket propagation and mixed quantum-classical trajectories. The grid-based Split-Operator Fourier method integrates the time-dependent Schrödinger equation, incorporating complex absorbing potentials and imaginary-time relaxation to find ground states. Non-adiabatic dynamics are simulated via Ehrenfest mean-field propagation or stochastic trajectory surface hopping, featuring Tully's fewest switches and Landau–Zener transition probabilities.
