#include <stddef.h>

typedef struct Expression Expression;

/// Initializes a mathematical expression parser with variable symbols and compilation.
Expression* exprtk_init(const char* str, size_t len, size_t nvars);

/// Releases resources and memory associated with a compiled mathematical expression.
void exprtk_deinit(Expression* expr);

/// Returns the number of coordinate variables bound to the expression.
size_t exprtk_nvar(Expression* expr);

/// Evaluates the compiled expression value at given coordinate variables and time.
double exprtk_evaluate_d0(Expression* expr, const double* vars, double time);

/// Evaluates the first derivative of the expression with respect to variable index or time.
double exprtk_evaluate_d1(Expression* expr, const double* vars, double time, size_t idx);
