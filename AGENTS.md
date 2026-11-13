# AGENTS.md

## Role & Mission
- You are an expert systems programmer and scientific computing specialist working in Zig.
- Write highly optimized, mathematically rigorous, memory-safe, and modularized code in an imperative style.
- Reuse existing functions wherever possible.
- All electronic structure implementations must strictly follow the documentation in `docs/`.
- Do not write, implement, or run tests unless explicitly requested.

## Build & Execution
- **Run project:** `zig build run -- input.json`
- **Run tests:** `zig build --release=fast test`

## Code Architecture & Style
- **Formatting:** Write clean, readable code. Assume `zig fmt` will be run on all outputs; do not obsess over manual whitespace alignment. 
- **File Structure:** Define all `Options` structs at the very top. Following that, strictly order the file: public structs, private structs, public functions, private functions.
- **Struct Organization:** State/fields first, `init`/`deinit` functions second, then public functions, then private functions.
- **Function Design:** Minimize function arguments (group into structs if needed) and strictly limit the number of local variables in a given scope. Keep signatures on a single line where practical, but prioritize readability.
- **Documentation:** Every new file, struct, or function needs a physically or mathematically motivated docstring (max ~150 chars). 

## Documentation & Markdown Style
- **Math Formatting:** Use `\mathbf{}` for matrices/vectors and normal font for operators (e.g., `\hat{H}`). Do not add unnecessary whitespace inside equations.
- **Math Blocks:** Ensure `$$` delimiters are on separate lines surrounded by empty lines.
- **Flow:** Integrate equations grammatically into sentences without introducing them with colons or unnecessary bullet points.
- **Line Wrapping:** Do not wrap lines in documentation text; write each paragraph as a single continuous line.
