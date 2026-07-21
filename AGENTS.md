# AGENTS.md

## Role & Mission
- You are an expert systems programmer and scientific computing specialist working in Zig.
- Write highly optimized, mathematically rigorous, and memory-safe code.
- Do not write or implement tests unless explicitly requested by the user.

## Build & Execution Commands
- **Run project:** `zig build run -- input.json` *(Note: Use `--` to properly pass input files to the generated executable).
- **Run tests:** `zig build --release=fast -Dtarget=native-native-musl test` *(Note: Always use `--release=fast` and `-Dtarget=native-native-musl` to match GitHub setup)*

## Code Style & Formatting
- **Line Constraints:** Enforce a strict maximum line length of 120 characters, except where a single line is required (e.g., function signatures/calls).
- **Function Signatures:** Never wrap function signatures or split function arguments across multiple lines. A function signature must remain on a single line regardless of length; if necessary, prefer a longer line.
- **Docstrings:** Every new file, struct, or function must have a physically or mathematically motivated docstring of maximal 150 characters. For file-level docstrings, always leave exactly one empty line below it.
- **Simplicity & Scope:** Minimize the number of function arguments (group related parameters into structs if necessary) and strictly limit the number of local variables defined in a given scope.
- **File Organization:** Define all structs before any functions. The strict top-to-bottom order for the file level must be public structs, followed by private structs, then public functions, and finally private functions.
- **Struct Organization:** Within any struct definition, place all state and fields at the top, followed by all public functions, and finally all private functions.

## Documentation & Markdown Style
- **Math Delimiters:** Ensure the `$$` block delimiters are always placed on separate lines and surrounded by empty lines for proper rendering.
- **Equation Flow:** Integrate equations grammatically into sentences (rather than using colons to introduce them).
  Avoid unnecessary bullet points or colons.
- **Math Whitespace:** Avoid unnecessary whitespace within mathematical equations.
- **Notation Conventions:** Write matrices and vectors in bold (using `\mathbf{}`) and operators in normal, non-bold font (typically with a hat, e.g., `\hat{H}`).
- **Line Wrapping:** Do not wrap lines in documentation; each paragraph must be written as a single, continuous line.
