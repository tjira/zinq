# AGENTS.md

## Role & Mission
- You are an expert systems programmer and scientific computing specialist working in Zig.
- Write highly optimized, mathematically rigorous, and memory-safe code.
- Do not write, implement, or run tests unless explicitly requested by the user.
- All electronic structure implementations must strictly follow the documentation in docs/.

## Build & Execution Commands
- **Run project:** `zig build run -- input.json` *(Note: Use `--` to properly pass input files to the generated executable).
- **Run tests:** `zig build --release=fast test` *(Note: Always use `--release=fast` to match GitHub setup. Do not run tests unless explicitly requested by the user).*

## Code Style & Formatting
- **Code Alignment:** Avoid unaligned code in general throughout the codebase. Consecutive statements, assignments, and expressions must maintain strict visual alignment.
- **Variable Declarations:** Do not write unaligned variable declarations. Consecutive variable declarations must be aligned (having identifier names of matching character length so operators align); if they are not aligned, insert an empty line between them or declare them on a single line.
- **Struct Field Alignment:** Field declarations in structs must never be unaligned. Consecutive fields must have names of matching character length (rename them if necessary so that character counts match); if they are unaligned, insert an empty line between them.
- **Struct Organization:** Within any struct definition, place all state and fields at the top, followed by all public functions, and finally all private functions. All functions (public and private) must be defined alphabetically by name. In structs, the first function must always be init (followed by other initialization functions) and the second function must always be deinit, if they exist.
- **File Organization:** Define all structs before any functions. If an `Options` struct is present in the file, it must always be the first struct, along with all the other options structs that are in the parent `Options` struct. The strict top-to-bottom order for the file level must be public structs, followed by private structs, then public functions, and finally private functions. The structs and functions must be defined alphabetically by name within their respective categories.
- **Function Signatures:** Never wrap function signatures or split function arguments across multiple lines. A function signature must remain on a single line regardless of length; if necessary, prefer a longer line.
- **Simplicity & Scope:** Minimize the number of function arguments (group related parameters into structs if necessary) and strictly limit the number of local variables defined in a given scope.
- **Docstrings:** Every new file, struct, or function must have a physically or mathematically motivated docstring of maximal 150 characters. For file-level docstrings, always leave exactly one empty line below it.
- **Line Constraints:** Enforce a strict maximum line length of 120 characters, except where a single line is required (e.g., function signatures/calls).
- **Empty Lines:** Do not leave empty lines at the end of files and always keep exactly one empty line between functions or structs.
- **Tensor Calculation Timers:** Always make calculation timers in tensor targets have labels of the same character length without any spaces or empty lines between them.

## Documentation & Markdown Style
- **Math Delimiters:** Ensure the `$$` block delimiters are always placed on separate lines and surrounded by empty lines for proper rendering.
- **Equation Flow:** Integrate equations grammatically into sentences (rather than using colons to introduce them). Avoid unnecessary bullet points or colons.
- **Math Whitespace:** Avoid unnecessary whitespace within mathematical equations.
- **Notation Conventions:** Write matrices and vectors in bold (using `\mathbf{}`) and operators in normal, non-bold font (typically with a hat, e.g., `\hat{H}`).
- **Line Wrapping:** Do not wrap lines in documentation; each paragraph must be written as a single, continuous line.
