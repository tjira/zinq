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
    <a href="https://github.com/tjira/zinq/releases/latest">
        <img src="https://img.shields.io/github/v/release/tjira/zinq?display_name=tag&style=for-the-badge"/>
    </a>
    <br>
    <a href="https://github.com/tjira/zinq">
        <img src="https://img.shields.io/github/languages/code-size/tjira/zinq?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/zinq">
        <img src="https://img.shields.io/endpoint?url=https://ghloc.vercel.app/api/tjira/zinq/badge?filter=.zig$&style=for-the-badge&format=human"/>
    </a>
    <a href="https://github.com/tjira/zinq/stargazers">
        <img src="https://img.shields.io/github/stars/tjira/zinq?style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/zinq/releases/latest">
        <img src="https://img.shields.io/github/downloads/tjira/zinq/total?style=for-the-badge"/>
    </a>
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
  Perform restricted or generalized Hartree-Fock calculation with the option for direct SCF.

* **Electronic Structure Analysis**  
  Compute energy derivatives and harmonic vibrational frequencies across supported methods.

### Time-Dependent Quantum Mechanics

* **Quantum Dynamics**  
  Simulate wavepacket dynamics in arbitrary dimensions and across multiple electronic states.

* **Surface Hopping**  
  Run nonadiabatic dynamics with Fewest Switches (FSSH) or Landau-Zener (LZSH) algorithms.

## Getting Zinq

### Prebuilt Releases

You can download the latest binaries from the [releases](https://github.com/tjira/zinq/releases/latest) page. The releases are provided for Linux, Windows and MacOS with the common CPU architectures. All binaries are statically linked with no external runtime dependencies. For less common platforms, see the [compilation](#Compilation) section.

### Compilation

Compiling Zinq is simple, running `make` will automatically download the Zig compiler to the project root and compile the Zinq binaries. The resulting executables are placed in the `zig-out` directory, organized by operating system and architecture. On Linux and Windows, most users will want the `x86_64` binary, while on MacOS the `aarch64` binary is usually appropriate. To verify the build, execute

```bash
./zig-out/<arch-os>/zinq
```

and check that the missing input message is displayed. If the message appears, the program is compiled correctly.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
