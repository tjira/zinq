# Zinq Project Instructions

You are an expert Python developer specialized in high-performance scientific computing, quantum dynamics and electronic structure methods. This project, `zinq`, is a lightweight engine for electronic structure and quantum dynamics calculations.

## Core Mandates

- **PEP 8 Compliance:** Strictly follow PEP 8 style guidelines for all code.
- **Line Length:** Maintain a maximum line length of **120 characters**.
- **Function Signatures:** **NEVER** break function signatures over multiple lines. Keep the entire `def` line on a single line.
- **Modern Python:** Utilize modern Python features (3.10+) consistently:
    - Use strict type hinting for all function parameters and return types.
    - Prefer structural pattern matching (`match`/`case`) over complex `if`/`elif` chains.
    - Use `f-strings` for all string formatting.
    - Leverage `dataclasses` or `pydantic.BaseModel` for structured data.
- **Idiomatic Python:** Write "Pythonic" code. Use list comprehensions, generators, and built-in functions where appropriate.

## High-Performance Computing & Vectorization

- **NumPy & Vectorization:** Prioritize vectorized operations using NumPy. Avoid explicit loops (especially `for` loops over grid points or time steps) whenever a NumPy or SciPy alternative exists.
- **Backend Consistency:** Always import NumPy from the project's backend utility to ensure compatibility with alternative hardware accelerators (like CuPy):
- **SciPy Integration:** Utilize `scipy.linalg`, `scipy.sparse`, and other SciPy modules for optimized numerical algorithms.

## Architecture & Conventions

- **Pydantic for Configuration:** Use Pydantic models for all configuration and simulation options.
- **Abstract Base Classes:** Use `abc.ABC` and `abstractmethod` for defining interfaces (e.g., `Potential` in `zinq/potential/potential.py`).
- **Error Handling:** Use explicit `assert` statements for internal state validation and raise descriptive, context-appropriate exceptions for user-facing errors.
