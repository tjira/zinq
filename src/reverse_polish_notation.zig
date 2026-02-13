//! Reverse Polish Notation (RPN) evaluation module.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");

const throw = error_handling.throw;

const STR2F = global_variables.STR2F;
const C2V = global_variables.C2V;
const RPN_MAX_STACK_SIZE = global_variables.RPN_MAX_STACK_SIZE;

/// Associativity types.
pub const Associativity = enum {
    Left,
    Right
};

/// Available operators.
pub const Operator = enum {
    Affirm,
    Plus,
    Minus,
    Negate,
    Multiply,
    Divide,
    Power
};

/// Reverse Polish Notation (RPN) structure.
pub fn ReversePolishNotation(comptime T: type) type {
    return struct {
        const Element = union(enum) {
            number: T,
            func: *const fn(T) T,
            op: Operator,
            variable: []const u8,
            constant: []const u8
        };

        data: std.ArrayList(Element),
        stack: std.ArrayList(T),

        /// Initializes the RPN structure.
        pub fn init(allocator: std.mem.Allocator) !@This() {
            return @This(){
                .data = std.ArrayList(Element){},
                .stack = try std.ArrayList(T).initCapacity(allocator, RPN_MAX_STACK_SIZE)
            };
        }

        /// Free the resources used by the RPN structure.
        pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
            self.data.deinit(allocator);
            self.stack.deinit(allocator);
        }

        /// Appends a number to the RPN output.
        pub fn append(self: *@This(), element: anytype, allocator: std.mem.Allocator) !void {
            if (@TypeOf(element) == T) {
                try self.data.append(allocator, Element{.number = element});
            } else if (@TypeOf(element) == *const fn(T) T) {
                try self.data.append(allocator, Element{.func = element});
            } else if (@TypeOf(element) == Operator) {
                try self.data.append(allocator, Element{.op = element});
            } else if (@TypeOf(element) == []const u8) {

                for (C2V.keys()) |constant| {
                    if (std.mem.eql(u8, element, constant)) {
                        try self.data.append(allocator, Element{.constant = element}); return;
                    }
                }

                try self.data.append(allocator, Element{.variable = element});
            } else {
                return throw(void, "INVALID ELEMENT TYPE FOR RPN", .{});
            }
        }

        /// Evaluates the RPN expression and returns the result.
        pub fn evaluate(self: *@This(), map: std.StringHashMap(T)) !T {
            for (self.data.items) |element| {
                switch (element) {
                    .number => |num| self.stack.appendAssumeCapacity(num),
                    .func => |func| self.stack.appendAssumeCapacity(func(self.stack.pop().?)),
                    .op => |op| switch (operatorArity(op)) {
                        1 => self.stack.appendAssumeCapacity(applyUnaryOperator(op, self.stack.pop().?)),
                        2 => self.stack.appendAssumeCapacity(applyBinaryOperator(op, self.stack.pop().?, self.stack.pop().?)),
                        else => unreachable
                    },
                    .variable => |variable| self.stack.appendAssumeCapacity(map.get(variable) orelse return throw(T, "UNKNOWN '{s}' VARIABLE IN RPN EVALUATION", .{variable})),
                    .constant => |constant| self.stack.appendAssumeCapacity(C2V.get(constant) orelse return throw(T, "UNKNOWN '{s}' CONSTANT IN RPN EVALUATION", .{constant})),
                }
            }

            return self.stack.pop().?;
        }

        /// Rerutns the RPN as a string.
        pub fn toString(self: @This(), allocator: std.mem.Allocator) ![]u8 {
            var buffer = std.ArrayList(u8){};

            outer: for (self.data.items, 0..) |element, i| {

                if (i > 0 and i < self.data.items.len) try buffer.print(allocator, " ", .{});

                switch (element) {

                    .number => |num| try buffer.print(allocator, "{d}", .{num}),

                    .func => |func| {

                        for (STR2F.keys(), STR2F.values()) |key, value| if (func == value) {
                            for (key) |c| try buffer.print(allocator, "{c}", .{std.ascii.toUpper(c)}); continue :outer;
                        };

                        return throw([]u8, "UNKNOWN FUNCTION IN RPN TO STRING", .{});
                    },

                    .op => |op| try buffer.print(allocator, "{s}", .{switch (op) {
                        .Plus => "+",
                        .Minus => "-",
                        .Affirm => "AFF",
                        .Negate => "NEG",
                        .Multiply => "*",
                        .Divide => "/",
                        .Power => "^",
                    }}),
                    
                    .variable => |variable| try buffer.print(allocator, "{s}", .{variable}),
                    .constant => |constant| try buffer.print(allocator, "{s}", .{constant}),
                }
            }

            return buffer.toOwnedSlice(allocator);
        }
    };
}

/// Applies the given unary operator to the provided number and returns the result.
pub fn applyUnaryOperator(op: Operator, n: anytype) @TypeOf(n) {
    return switch (op) {
        .Affirm => n,
        .Negate => -n,
        else => unreachable
    };
}

/// Applies the given binary operator to the two provided numbers and returns the result.
pub fn applyBinaryOperator(op: Operator, n1: anytype, n2: @TypeOf(n1)) @TypeOf(n1, n2) {
    return switch (op) {
        .Plus => n2 + n1,
        .Minus => n2 - n1,
        .Multiply => n2 * n1,
        .Divide => n2 / n1,
        .Power => std.math.pow(@TypeOf(n1, n2), n2, n1),
        else => unreachable
    };
}

/// Returns te associativity of the given operator.
pub fn operatorAssociativity(op: Operator) Associativity {
    return switch (op) {
        .Plus, .Minus, .Multiply, .Divide => .Left,
        .Power, .Negate, .Affirm => .Right
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

/// Returns the operator arity (number of operands) for the given operator.
pub fn operatorArity(op: Operator) u8 {
    return switch (op) {
        .Plus, .Minus, .Multiply, .Divide, .Power => 2,
        .Negate, .Affirm => 1,
    };
}

/// Returns the precedence of the given operator.
pub fn operatorPrecedence(op: Operator) u8 {
    return switch (op) {
        .Plus, .Minus => 2,
        .Multiply, .Divide => 3,
        .Power, .Negate, .Affirm => 4,
    };
}
