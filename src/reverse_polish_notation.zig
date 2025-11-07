//! Reverse Polish Notation (RPN) evaluation module.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");

const throw = error_handling.throw;

const RPN_MAX_STACK_SIZE = global_variables.RPN_MAX_STACK_SIZE;

/// Associativity types.
pub const Associativity = enum {
    Left,
    Right
};

/// Available operators.
pub const Operator = enum {
    Plus,
    Minus,
    Multiply,
    Divide,
    Power
};

/// Reverse Polish Notation (RPN) structure.
pub fn ReversePolishNotation(comptime T: type) type {
    return struct {
        const Element = union(enum) {
            number: T,
            op: Operator,
            variable: []const u8
        };

        data: std.ArrayList(Element),
        stack: std.ArrayList(T),
        allocator: std.mem.Allocator,

        /// Initializes the RPN structure.
        pub fn init(allocator: std.mem.Allocator) !@This() {
            return @This(){
                .data = std.ArrayList(Element){},
                .stack = try std.ArrayList(T).initCapacity(allocator, RPN_MAX_STACK_SIZE),
                .allocator = allocator
            };
        }

        /// Free the resources used by the RPN structure.
        pub fn deinit(self: *@This()) void {
            self.data.deinit(self.allocator);
            self.stack.deinit(self.allocator);
        }

        /// Appends a number to the RPN output.
        pub fn append(self: *@This(), element: anytype) !void {
            if (@TypeOf(element) == T) {
                try self.data.append(self.allocator, Element{.number = element});
            } else if (@TypeOf(element) == Operator) {
                try self.data.append(self.allocator, Element{.op = element});
            } else if (@TypeOf(element) == []const u8) {
                try self.data.append(self.allocator, Element{.variable = element});
            } else {
                return throw(void, "INVALID ELEMENT TYPE FOR RPN", .{});
            }
        }

        /// Evaluates the RPN expression and returns the result.
        pub fn evaluate(self: *@This(), map: std.StringHashMap(T)) !T {
            for (self.data.items) |element| {
                switch (element) {
                    .number => |num| self.stack.appendAssumeCapacity(num),
                    .op => |op| self.stack.appendAssumeCapacity(applyOperator(op, self.stack.pop().?, self.stack.pop().?)),
                    .variable => |variable| self.stack.appendAssumeCapacity(map.get(variable) orelse return throw(T, "UNKNOWN '{s}' VARIABLE IN RPN EVALUATION", .{variable})),
                }
            }

            return self.stack.pop().?;
        }
    };
}

/// Applys the given operator to two operands.
pub fn applyOperator(op: Operator, n1: anytype, n2: @TypeOf(n1)) @TypeOf(n1, n2) {
    return switch (op) {
        .Plus => n2 + n1,
        .Minus => n2 - n1,
        .Multiply => n2 * n1,
        .Divide => n2 / n1,
        .Power => std.math.pow(@TypeOf(n1, n2), n2, n1)
    };
}

/// Returns te associativity of the given operator.
pub fn operatorAssociativity(op: Operator) Associativity {
    return switch (op) {
        .Plus, .Minus, .Multiply, .Divide => .Left,
        .Power => .Right
    };
}

/// Returns the operator corresponding to the given character.
pub fn operatorFromChar(char: u8) !Operator {
    return switch (char) {
        '+' => Operator.Plus,
        '-' => Operator.Minus,
        '*' => Operator.Multiply,
        '/' => Operator.Divide,
        '^' => Operator.Power,
        else => return throw(Operator, "UNKNOWN OPERATOR CHARACTER", .{}),
    };
}

/// Returns the precedence of the given operator.
pub fn operatorPrecedence(op: Operator) u8 {
    return switch (op) {
        .Plus, .Minus => 2,
        .Multiply, .Divide => 3,
        .Power => 4
    };
}
