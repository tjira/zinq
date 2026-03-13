//! File with custom electronic potential implementation using the Shunting Yard algorithm.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const reverse_polish_notation = @import("reverse_polish_notation.zig");
const shunting_yard = @import("shunting_yard.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const ReversePolishNotation = reverse_polish_notation.ReversePolishNotation;

const shuntingYard = shunting_yard.shuntingYard;

/// Struct holding parameters for the custom electronic potential that uses the Shunting Yard algorithm to parse mathematical expressions.
pub fn CustomPotential(comptime T: type) type {
    return struct {
        variables: []const []const u8,
        matrix: []const []const []const u8,

        data: ?struct {
            rpn_array: []ReversePolishNotation(T),
        } = null,

        /// Diabatic potential evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), time: T) !void {
            for (0..U.rows) |i| for (i..U.cols) |j| {
                U.ptr(i, j).* = try evaluateDiabaticElement(self, i, j, position, time); U.ptr(j, i).* = U.at(i, j);
            };
        }

        /// Diabatic potential matrix element evaluator.
        pub fn evaluateDiabaticElement(self: @This(), i: usize, j: usize, position: RealVector(T), time: T) !T {
            if (i >= self.matrix.len or j >= self.matrix.len) {

                std.log.err("INVALID INDEX ({d}, {d}) WHEN EVALUATING DIABATIC MATRIX ELEMENT, THE POTENTIAL MATRIX IS {d}X{d}", .{i, j, 2, 2});

                return error.ProgrammingError;
            }

            var buffer: [1024]u8 = undefined; var fba = std.heap.FixedBufferAllocator.init(&buffer);

            var map = std.StringHashMap(T).init(fba.allocator()); defer map.deinit();

            try map.put("t", time);

            for (0..self.variables.len) |k| try map.put(self.variables[k], position.at(k));

            if (self.data) |cpd| return try cpd.rpn_array[i * self.matrix.len + j].evaluate(map);

            std.log.err("CUSTOM POTENTIAL DATA NOT INITIALIZED, CALL init() BEFORE EVALUATING THE POTENTIAL", .{});

            return error.ProgrammingError;
        }

        /// Parse and initialize the custom potential from the provided expression matrix.
        pub fn init(self: *@This(), allocator: std.mem.Allocator) !void {
            var rpn_array = try allocator.alloc(ReversePolishNotation(T), self.matrix.len * self.matrix.len);

            for (0..self.matrix.len) |i| if (self.matrix[i].len != self.matrix.len) {

                std.log.err("THE POTENTIAL MATRIX MUST BE SQUARE, BUT THE PROVIDED MATRIX HAS {d} ROWS AND {d} COLUMNS", .{self.matrix.len, self.matrix[i].len});

                return error.InputError;
            };

            for (0..self.matrix.len) |i| for (0..self.matrix.len) |j| {
                rpn_array[i * self.matrix.len + j] = try shuntingYard(T, self.matrix[i][j], self.variables, allocator);
            };

            self.data = .{
                .rpn_array = rpn_array,
            };
        }

        pub fn deinit(self: @This(), allocator: std.mem.Allocator) void {
            if (self.data) |data| {

                for (data.rpn_array) |*rpn| rpn.deinit(allocator);

                allocator.free(data.rpn_array);
            }
        }

        /// Custom potential JSON parser.
        pub fn jsonParseFromValue(allocator: std.mem.Allocator, source: std.json.Value, options: std.json.ParseOptions) !@This() {
            if (source != .object) return error.UnexpectedToken;

            const variables_val = source.object.get("variables") orelse return error.MissingField;
            const matrix_val = source.object.get("matrix") orelse return error.MissingField;

            const parsed_variables = try std.json.innerParseFromValue([]const []const u8, allocator, variables_val, options);
            const parsed_matrix = try std.json.innerParseFromValue([]const []const []const u8, allocator, matrix_val, options);

            return .{.variables = parsed_variables, .matrix = parsed_matrix};
        }

        /// Custom potential JSON stringifier.
        pub fn jsonStringify(self: *const @This(), jws: anytype) !void {
            try jws.beginObject();
            
            try jws.objectField("variables");
            try jws.write(self.variables);
            
            try jws.objectField("matrix");
            try jws.write(self.matrix);
            
            try jws.endObject();
        }
    };
}
