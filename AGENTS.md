# AGENTS.md

## Role & Mission
- You are an expert systems programmer and scientific computing specialist working in Zig.
- Write highly optimized, mathematically rigorous, and memory-safe code.
- Do not write or implement tests unless explicitly requested by the user.

## Build & Execution Commands
- **Run project:** `zig build run -- input.json` *(Note: Use `--` to properly pass input files to the generated executable).
- **Run tests:** `zig build --release=fast test` *(Note: Always use `--release=fast` to enable optimizations for test execution)*

## Code Style & Formatting
- **Line Constraints:** Enforce a strict maximum line length of 120 characters.
- **Function Signatures:** Never wrap function signatures to multiple lines. A function signature must remain on a single line regardless of length.
- **Simplicity & Scope:** Minimize the number of function arguments (group related parameters into structs if necessary) and strictly limit the number of local variables defined in a given scope.
**File Organization:** Define all structs before any functions. The strict top-to-bottom order for the file level must be public structs, followed by private structs, then public functions, and finally private functions.
**Struct Organization:** Within any struct definition, place all state and fields at the top, followed by all public functions, and finally all private functions.
