//! Reverse Polish Notation (RPN) evaluation module.

const std = @import("std");

const global_variables = @import("global_variables.zig");

const STR2F = global_variables.STR2F;
const C2V = global_variables.C2V;

/// Associativity types.
pub const Associativity = enum { Left, Right };

/// Available operators.
pub const Operator = enum { Affirm, Plus, Minus, Negate, Multiply, Divide, Power };

/// Reverse Polish Notation (RPN) structure.
pub fn ReversePolishNotation(comptime T: type) type {
    return struct {
        const Element = union(enum) { number: T, func: *const fn (T) T, op: Operator, variable: []const u8, constant: []const u8 };

        data: std.ArrayList(Element),

        /// Initializes the RPN structure.
        pub fn init() !@This() {
            return @This(){
                .data = std.ArrayList(Element){},
            };
        }

        /// Free the resources used by the RPN structure.
        pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
            self.data.deinit(allocator);
        }

        /// Appends a number to the RPN output.
        pub fn append(self: *@This(), element: anytype, allocator: std.mem.Allocator) !void {
            if (@TypeOf(element) == T) {
                try self.data.append(allocator, Element{ .number = element });
            } else if (@TypeOf(element) == *const fn (T) T) {
                try self.data.append(allocator, Element{ .func = element });
            } else if (@TypeOf(element) == Operator) {
                try self.data.append(allocator, Element{ .op = element });
            } else if (@TypeOf(element) == []const u8) {
                for (C2V.keys()) |constant| {
                    if (std.mem.eql(u8, element, constant)) {
                        try self.data.append(allocator, Element{ .constant = element });
                        return;
                    }
                }

                try self.data.append(allocator, Element{ .variable = element });
            } else {
                std.log.err("UNSUPPORTED ELEMENT TYPE '{s}' IN RPN APPEND, EXPECTED A NUMBER, FUNCTION POINTER, OPERATOR OR VARIABLE/CONSTANT NAME", .{std.meta.typeName(@TypeOf(element))});

                return error.InvalidInput;
            }
        }

        /// Evaluates the RPN expression and returns the result.
        pub fn evaluate(self: @This(), map: std.StringHashMap(T)) !T {
            var stack: [2048]T = undefined;
            var length: usize = 0;

            for (self.data.items) |element| {
                switch (element) {
                    .number => |num| {
                        if (length >= stack.len) {
                            std.log.err("RPN STACK OVERFLOW, MAX STACK SIZE IS {d}, MAXIMUM STACK SIZE CAN BE CHANGED IN THE SOURCE CODE", .{stack.len});

                            return error.ProgrammingError;
                        }
                        stack[length] = num;
                        length += 1;
                    },
                    .func => |func| {
                        if (length < 1) {
                            std.log.err("RPN STACK UNDERFLOW WHEN APPLYING FUNCTION, NOT ENOUGH OPERANDS IN THE STACK", .{});

                            return error.ProgrammingError;
                        }
                        stack[length - 1] = func(stack[length - 1]);
                    },
                    .op => |op| switch (operatorArity(op)) {
                        1 => {
                            if (length < 1) {
                                std.log.err("RPN STACK UNDERFLOW WHEN APPLYING UNARY OPERATOR, NOT ENOUGH OPERANDS IN THE STACK", .{});

                                return error.ProgrammingError;
                            }
                            stack[length - 1] = applyUnaryOperator(op, stack[length - 1]);
                        },
                        2 => {
                            if (length < 2) {
                                std.log.err("RPN STACK UNDERFLOW WHEN APPLYING BINARY OPERATOR, NOT ENOUGH OPERANDS IN THE STACK", .{});

                                return error.ProgrammingError;
                            }
                            const n1 = stack[length - 1];
                            const n2 = stack[length - 2];
                            stack[length - 2] = applyBinaryOperator(op, n1, n2);
                            length -= 1;
                        },
                        else => {
                            std.log.err("UNKNOWN OPERATOR ARITY {d} FOR OPERATOR", .{operatorArity(op)});

                            return error.ProgrammingError;
                        },
                    },
                    .variable => |variable| {
                        if (map.get(variable) == null) {
                            std.log.err("UNDEFINED VARIABLE '{s}' IN RPN EVALUATION", .{variable});

                            return error.InvalidInput;
                        }

                        if (length >= stack.len) {
                            std.log.err("RPN STACK OVERFLOW, MAX STACK SIZE IS {d}", .{stack.len});

                            return error.ProgrammingError;
                        }

                        stack[length] = map.get(variable).?;
                        length += 1;
                    },
                    .constant => |constant| {
                        if (C2V.get(constant) == null) {
                            std.log.err("UNDEFINED CONSTANT '{s}' IN RPN EVALUATION", .{constant});

                            return error.InvalidInput;
                        }

                        if (length >= stack.len) {
                            std.log.err("RPN STACK OVERFLOW, MAX STACK SIZE IS {d}", .{stack.len});

                            return error.ProgrammingError;
                        }

                        stack[length] = C2V.get(constant).?;
                        length += 1;
                    },
                }
            }

            if (length != 1) {
                std.log.err("ERROR IN EVALUATION OF RPN EXPRESSION, EXPECTED A SINGLE RESULT BUT GOT {d} ELEMENTS IN THE STACK", .{length});

                return error.ProgrammingError;
            }

            return stack[0];
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
                            for (key) |c| try buffer.print(allocator, "{c}", .{std.ascii.toUpper(c)});
                            continue :outer;
                        };

                        std.log.err("UNKNOWN FUNCTION PASSED TO RPN TO STRING", .{});

                        return error.InvalidInput;
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
        else => unreachable,
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
        else => unreachable,
    };
}

/// Returns te associativity of the given operator.
pub fn operatorAssociativity(op: Operator) Associativity {
    return switch (op) {
        .Plus, .Minus, .Multiply, .Divide => .Left,
        .Power, .Negate, .Affirm => .Right,
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
        else => return error.InvalidInput,
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
