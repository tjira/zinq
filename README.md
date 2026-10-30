<h1 align="center">Zinq</h1>

<h4 align="center">
  <a href="#Features">Features</a>
  ·
  <a href="#Citation">Citation</a>
  ·
  <a href="https://tjira.github.io/zinq">Docs</a>
</h4>

<p align="center">
    <a href="https://github.com/tjira/zinq/pulse">
        <img src="https://img.shields.io/github/last-commit/tjira/zinq?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/zinq/blob/master/LICENSE">
        <img src="https://img.shields.io/github/license/tjira/zinq?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/zinq/actions/workflows/test.yml">
        <img src="https://img.shields.io/github/actions/workflow/status/tjira/zinq/test.yml?style=for-the-badge&label=test"/>
    </a>
    <a href="https://pypi.org/project/zinq">
        <img src="https://img.shields.io/pypi/v/zinq?style=for-the-badge"/>
    </a>
    <br>
    <a href="https://github.com/tjira/zinq">
        <img src="https://img.shields.io/github/languages/code-size/tjira/zinq?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/zinq">
        <img src="https://img.shields.io/endpoint?url=https://ghloc.vercel.app/api/tjira/zinq/badge?filter=.zig&style=for-the-badge&format=human"/>
    </a>
    <a href="https://github.com/tjira/zinq/stargazers">
        <img src="https://img.shields.io/github/stars/tjira/zinq?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/zinq/releases/latest">
        <img src="https://img.shields.io/github/downloads/tjira/zinq/total?style=for-the-badge"/>
    </a>
    <br>
</p>

<p align="center">
Lightweight Zig framework for electronic structure theory, quantum chemistry, and mathematical algorithms. Written from scratch, it favors simple design while relying on efficient algorithms.
</p>

## Features

Zinq provides tools for both time-independent and time-dependent quantum mechanical simulations.

### Time-Independent Quantum Mechanics

* **Integrals over Gaussian Basis Functions**  
  Compute 1- and 2-electron integrals and geometric derivative integrals.

* **Hartree–Fock Methods**  
  Perform restricted or generalized Hartree–Fock calculations with DIIS acceleration.

* **Density Functional Theory**  
  Evaluate Kohn–Sham DFT using LDA, GGA, meta-GGA, and hybrid functionals on molecular grids.

* **Post-Hartree–Fock Methods**  
  Perform arbitrary-order Møller–Plesset perturbation theory and arbitrary-excitation CI.

* **Multireference Methods**  
  Perform CASSCF calculations with simultaneous orbital and CI optimization.

* **Nuclear Derivatives & Vibrational Analysis**  
  Compute analytic or numerical gradients, numerical Hessians, and harmonic vibrational frequencies.

* **Molecular Geometry Optimization**  
  Optimize molecular structures using BFGS quasi-Newton and steepest descent algorithms.

* **Population Analysis & Bond Orders**  
  Calculate Mulliken and Löwdin populations, and Mayer and Wiberg bond orders.

### Time-Dependent Quantum Mechanics

* **Quantum Wavepacket Dynamics**  
  Simulate multidimensional wavepacket dynamics using split-operator FFT propagation.

* **Trajectory Surface Hopping**  
  Run nonadiabatic dynamics using fewest switches (FSSH) and Landau–Zener algorithms.

* **Ehrenfest Dynamics**  
  Perform mean-field nonadiabatic simulations with continuous electronic-nuclear coupling.

* **Classical Molecular Dynamics**  
  Propagate trajectories with thermostats and sample Maxwell–Boltzmann or Wigner ensembles.

## Citation

If you use Zinq in your research, please cite the project on [Zenodo](https://doi.org/10.5281/zenodo.18386143) using the following general BibTeX entry. If you are referring to a specific version, please visit the Zenodo page and download the corresponding citation there.

```bibtex
@software{ZinqJira2026,
    author    = {Tomáš Jíra},
    title     = {tjira/zinq},
    year      = 2026,
    publisher = {Zenodo},
    doi       = {10.5281/zenodo.18386143},
    url       = {https://doi.org/10.5281/zenodo.18386143},
}
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

<p align="left">
    <a href="https://github.com/tjira/zinq">
        <img src="https://img.shields.io/badge/main-github-black?logo=github&style=for-the-badge"/>
    </a>
    <a href="https://gitea.com/tjira/zinq">
        <img src="https://img.shields.io/badge/mirror-gitea-green?logo=gitea&style=for-the-badge"/>
    </a>
    <a href="https://gitlab.com/tojira/zinq">
        <img src="https://img.shields.io/badge/mirror-gitlab-orange?logo=gitlab&style=for-the-badge"/>
    </a>
    <a href="https://codeberg.org/tjira/zinq">
        <img src="https://img.shields.io/badge/mirror-codeberg-blue?logo=codeberg&style=for-the-badge"/>
    </a>
</p>
