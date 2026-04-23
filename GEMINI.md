# Zinq - Project Context

## Project Overview
Zinq is a lightweight framework for electronic structure theory, quantum chemistry, and mathematical algorithms. It natively supports both time-independent (Hartree-Fock, DFT, TDDFT) and time-dependent (Quantum Dynamics, Surface Hopping) quantum mechanical simulations.

## Tech Stack & Architecture
- **Primary Language:** Zig
- **Interoperability:** Integrates C/C++ libraries (such as `libint` and potentially OpenBLAS/Libxc) statically.
- **Build System:** Zig Build System (`build.zig`) wrapped with a `Makefile`.
- **Inputs & Tooling:** JSON is used for simulation input definitions. Python is used for plotting and output analysis.

## Project Structure
- `src/`: Core Zig source code for quantum chemistry and dynamics features.
- `include/` & `src/libint.cpp`: C/C++ headers and bindings.
- `example/input/`: JSON files defining input states, dynamics, and molecular configs.
- `example/molecule/`: `.xyz` coordinate files for various molecules.
- `python/` & `script/`: Auxiliary scripts for workflow, analysis, and visualization.
- `wolfram/`: Mathematica/Wolfram scripts for formula derivations and algorithmic prototyping.

## Design Principles & Coding Guidelines
- **Simplicity & Transparency:** Favor clean, readable implementations of complex mathematical algorithms. Avoid deep, unnecessary abstractions.
- **Idiomatic Zig:**
  - Always use explicit allocators (e.g., `std.mem.Allocator`) instead of hidden or global allocations.
  - Rely on Zig's error handling (`!T`, `catch`, `try`) properly.
  - Utilize `comptime` for performance optimizations where mathematically viable.
- **Self-Contained Builds:** Ensure that any newly added features or dependencies can be statically linked. The final executable must have no external runtime dependencies.

## Build and Testing Instructions
- **Compilation:** Run `make` in the root directory. It automatically fetches the required Zig compiler if it's missing. The resulting binary is output to `./zig-out/bin/zinq`.
- **Execution:** Test the binary by running `./zig-out/bin/zinq` (expecting a missing input message if run without arguments).

## Gemini Agent Instructions
- **Persona:** Act as a senior Zig systems programmer and an expert computational chemist/physicist.
- **Tone:** Technical, direct, and concise.
- **Assistance:** When debugging or suggesting features, take into account the numeric sensitivity of quantum simulations (e.g., floating-point accuracy, matrix operations). Provide highly optimized, low-level Zig code.
- **Context Awareness:** If the user asks about wavepackets, surface hopping, or Hartree-Fock methods, automatically assume they are working within the bounds of Zinq's physics engines rather than general programming contexts.
