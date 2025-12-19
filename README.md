<h1 align="center">Zinq</h1>

<h4 align="center">
  <a href="https://github.com/tjira/zinq#features">Features</a>
  ·
  <a href="https://github.com/tjira/zinq#compilation">Compilation</a>
  ·
  <a href="https://tjira.github.io/zinq/">Docs</a>
</h4>

<p align="center">
    <a href="https://github.com/tjira/zinq/pulse">
        <img src="https://img.shields.io/github/last-commit/tjira/zinq?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/zinq/blob/master/license">
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
        <img src="https://img.shields.io/endpoint?url=https://ghloc.vercel.app/api/tjira/zinq/badge?filter=.zig,!hermite_quadrature_nodes.zig&style=for-the-badge&format=human"/>
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
A lightweight Zig framework for electronic structure theory, quantum chemistry, and mathematical algorithms. Written from scratch, it favors simple design and transparent implementation while relying on efficient algorithms.
</p>

## Features

Zinq provides tools for both time-independent and time-dependent quantum mechanical simulations.

### Time-Independent Quantum Mechanics

* **Integrals over Gaussian Basis Functions**  
  Compute integrals over Gaussian basis functions from .xyz geometries and basis files.

* **Hartree–Fock Methods**  
  Perform restricted or generalized Hartree-Fock calculation with DIIS accelerator.

* **Post-Hartree–Fock Methods**  
  Use variety of selected perturbative or variational post-Hartree-Fock methods.

* **Electronic Structure Analysis**  
  Compute energy derivatives and harmonic vibrational frequencies across supported methods.

### Time-Dependent Quantum Mechanics

* **Quantum Dynamics**  
  Simulate wavepacket dynamics in arbitrary dimensions and across multiple electronic states.

* **Dirac–Frenkel Variational Principle**  
  Propagate a parametrized wavefunction using the Dirac–Frenkel variational principle.

* **Surface Hopping**  
  Run nonadiabatic dynamics with various surface hopping algorithms.

## Getting Zinq

### Prebuilt Releases

You can download the latest binaries from the [releases](https://github.com/tjira/zinq/releases/latest) page. The releases are provided for Linux, Windows and MacOS with the common CPU architectures. All binaries are statically linked with no external runtime dependencies. For less common platforms, see the [compilation](#Compilation) section. The binaries can also be installed using `pip` from [PyPI](https://pypi.org/project/zinq).

### Compilation

Compiling Zinq is simple, running `make` will automatically download the Zig compiler to the project root and compile the Zinq binaries. The resulting executables are placed in the `zig-out` directory, organized by operating system and architecture. On Linux and Windows, most users will want the `x86_64` binary, while on MacOS the `aarch64` binary is usually appropriate. To verify the build, execute

```bash
./zig-out/<arch-os>/zinq
```

and check that the missing input message is displayed. If the message appears, the program is compiled correctly.

## Python Wrappers

You can access Zinq's functionality through the dedicated Python wrappers. If you installed Zinq via [PyPI](https://pypi.org/project/zinq), these wrappers are pre-installed and ready to use. However, if you are compiling from source, you will need to build them manually. From the project root, simply run `make pip` to build the wheel and install the package (or `make wheel` to build without installing). Once set up, the `zinq` command is added to your PATH. You also gain access to:

* `hf` - Runs Hartree-Fock calculations on `.xyz` geometry files.
* `molint` - Computes molecular integrals from `.xyz` geometry files.
* `mp2` - Performs MP2 energy calculations on `.xyz` geometry files.
* `prime` - Generates prime numbers as a standalone utility tool.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

<p align="left">
    <a href="https://gitlab.com/tojira/zinq">
        <img src="https://img.shields.io/badge/mirror-gitlab-orange?logo=gitlab&style=for-the-badge"/>
    </a>
    <a href="https://codeberg.org/tjira/zinq">
        <img src="https://img.shields.io/badge/mirror-codeberg-blue?logo=codeberg&style=for-the-badge"/>
    </a>
</p>
