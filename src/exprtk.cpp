#include <exprtk.hpp>

struct Expression {
    exprtk::expression<double> expression;
    exprtk::symbol_table<double>  symbols;
    std::vector<double> vars; double time;
};

extern "C" {
    /// Initializes a mathematical expression parser with variable symbols and compilation.
    Expression* exprtk_init(const char* str, size_t len, size_t nvar) {
        Expression* expr = new Expression(); expr->vars.resize(nvar);

        for (size_t i = 1; i <= nvar; ++i) {
            expr->symbols.add_variable("r" + std::to_string(i), expr->vars.at(i - 1));
        }

        expr->symbols.add_variable("t", expr->time); expr->symbols.add_constants();

        expr->expression.register_symbol_table(expr->symbols);

        if (!exprtk::parser<double>().compile(std::string(str, len), expr->expression)) {
            delete expr; return nullptr;
        }

        return expr;
    }

    /// Releases resources and memory associated with a compiled mathematical expression.
    void exprtk_deinit(Expression* expr) {
        if (!expr) return; delete expr;
    }

    /// Evaluates the compiled expression value at given coordinate variables and time.
    double exprtk_evaluate_d0(Expression* expr, const double* vars, double time) {
        if (!expr) return 0;

        std::copy(vars, vars + expr->vars.size(), expr->vars.begin());

        expr->time = time;

        return expr->expression.value();
    }

    /// Evaluates the first derivative of the expression with respect to variable index or time.
    double exprtk_evaluate_d1(Expression* expr, const double* vars, double time, size_t idx) {
        if (!expr or idx > expr->vars.size()) return 0;

        std::copy(vars, vars + expr->vars.size(), expr->vars.begin());

        expr->time = time;

        if (idx == expr->vars.size()) {
            return exprtk::derivative(expr->expression, expr->time);
        }

        return exprtk::derivative(expr->expression, expr->vars[idx]);
    }

    /// Returns the number of coordinate variables bound to the expression.
    size_t exprtk_nvar(Expression* expr) {
        if (!expr) return 0; return expr->vars.size();
    }
}
