#include <stddef.h>

typedef struct Expression Expression;

Expression* exprtk_init(const char* str, size_t len, size_t nvars); void exprtk_deinit(Expression* expr);

size_t exprtk_nvar(Expression* expr);

double exprtk_evaluate_d0(Expression* expr, const double* vars, double time            );
double exprtk_evaluate_d1(Expression* expr, const double* vars, double time, size_t idx);
