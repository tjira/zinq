//! Implementation of the Shunting Yard algorithm in Zig.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const reverse_polish_notation = @import("reverse_polish_notation.zig");

const Operator = reverse_polish_notation.Operator;
const ReversePolishNotation = reverse_polish_notation.ReversePolishNotation;

const operatorAssociativity = reverse_polish_notation.operatorAssociativity;
const operatorFromChar = reverse_polish_notation.operatorFromChar;
const operatorPrecedence = reverse_polish_notation.operatorPrecedence;
const throw = error_handling.throw;

/// Parses the input expression using the Shunting Yard algorithm and returns the RPN output.
pub fn shuntingYard(comptime T: type, input: []const u8, variables: []const []const u8, allocator: std.mem.Allocator) !ReversePolishNotation(T) {
    var rpn = try ReversePolishNotation(T).init(allocator);

    var stack = std.ArrayList(union(enum) {op: Operator, bracket: u8}){}; defer stack.deinit(allocator);

    var i: usize = 0; var j: usize = 0; var buffer: [64]u8 = undefined;

    parser: while (i < input.len) : (i += 1) {

        if (std.ascii.isWhitespace(input[i])) continue;

        if (std.ascii.isDigit(input[i])) {

            while (i < input.len and (std.ascii.isDigit(input[i]) or input[i] == '.')) : (i += 1) {
                buffer[j] = input[i]; j += 1;
            }

            const number = try std.fmt.parseFloat(T, buffer[0..j]);

            try rpn.append(number, allocator); i -= 1; j = 0; continue;
        }

        else if (input[i] == 'e') {

            const number: T = std.math.e;

            try rpn.append(number, allocator); j = 0; continue;
        }

        else if (input[i] == '+' or input[i] == '-' or input[i] == '*' or input[i] == '/' or input[i] == '^') {

            const op = try operatorFromChar(input[i]); const opp = operatorPrecedence(op); const opa = operatorAssociativity(op);

            while (stack.items.len > 0 and switch(stack.getLast()) {.op => true, .bracket => |bracket| bracket != '('}) {
                if (operatorPrecedence(stack.getLast().op) > opp or (operatorPrecedence(stack.getLast().op) == opp and opa == .Left)) try rpn.append(stack.pop().?.op, allocator) else break;
            }

            try stack.append(allocator, .{.op = op});
        }

        else if (input[i] == '(') {
            try stack.append(allocator, .{.bracket = input[i]});
        }

        else if (input[i] == ')') {

            while (stack.items.len > 0 and switch(stack.getLast()) {.op => true, .bracket => |bracket| bracket != '('}) {
                try rpn.append(stack.pop().?.op, allocator);
            }

            if (stack.items.len == 0) return throw(ReversePolishNotation(T), "MISMATCHED PARENTHESES IN EXPRESSION", .{});

            _ = stack.pop();
        }

        else {

            for (variables) |variable| {
                if (std.mem.eql(u8, variable, input[i..i + variable.len])) {
                    try rpn.append(variable, allocator); i += variable.len - 1; continue :parser;
                }
            }

            return throw(ReversePolishNotation(T), "UNKNOWN TOKEN IN EXPRESSION", .{});
        }
    }

    while (stack.items.len > 0) try rpn.append(stack.pop().?.op, allocator);

    return rpn;
}
