//! File with custom electronic potential implementation using the Shunting Yard algorithm.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const reverse_polish_notation = @import("reverse_polish_notation.zig");
const shunting_yard = @import("shunting_yard.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const ReversePolishNotation = reverse_polish_notation.ReversePolishNotation;

const shuntingYard = shunting_yard.shuntingYard;
const throw = error_handling.throw;

var custom_potential_data: ?CustomPotentialData(f64) = null;

/// Global struct holding the data for the custom electronic potential.
pub fn CustomPotentialData(comptime T: type) type {
    return struct {
        rpn_array: []ReversePolishNotation(T),
        map: std.StringHashMap(T),
        allocator: std.mem.Allocator,

        /// Free the resources.
        pub fn deinit(self: *@This()) void {
            for (self.rpn_array) |*rpn| rpn.deinit();
            self.allocator.free(self.rpn_array);
            self.map.deinit();
        }
    };
}

/// Struct holding parameters for the custom electronic potential that uses the Shunting Yard algorithm to parse mathematical expressions.
pub fn CustomPotential(comptime T: type) type {
    return struct {
        variables: []const []const u8,
        matrix: []const []const []const u8,

        /// Diabatic potential evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), time: T) !void {
            custom_potential_data.?.map.putAssumeCapacity("t", time);

            for (0..position.len) |q| custom_potential_data.?.map.putAssumeCapacity(self.variables[q], position.at(q));

            for (0..U.rows) |i| for (i..U.cols) |j| {
                U.ptr(i, j).* = try custom_potential_data.?.rpn_array[i * U.cols + j].evaluate(custom_potential_data.?.map); U.ptr(j, i).* = U.at(i, j);
            };
        }

        /// Diabatic potential matrix element evaluator.
        pub fn evaluateDiabaticElement(self: @This(), i: usize, j: usize, position: RealVector(T), time: T) !T {
            if (i >= self.matrix.len or j >= self.matrix.len) return throw(T, "INVALID INDEX WHEN EVALUATING DIABATIC MATRIX ELEMENT", .{});

            custom_potential_data.?.map.putAssumeCapacity("t", time);

            for (0..position.len) |q| custom_potential_data.?.map.putAssumeCapacity(self.variables[q], position.at(q));

            return try custom_potential_data.?.rpn_array[i * self.matrix.len + j].evaluate(custom_potential_data.?.map);
        }

        /// Parse and initialize the custom potential from the provided expression matrix.
        pub fn init(self: @This(), allocator: std.mem.Allocator) !?CustomPotentialData(T) {
            var rpn_array = try allocator.alloc(ReversePolishNotation(T), self.matrix.len * self.matrix.len);
            var map = std.StringHashMap(T).init(allocator);

            for (0..self.matrix.len) |i| if (self.matrix[i].len != self.matrix.len) return throw(?CustomPotentialData(T), "CUSTOM POTENTIAL MATRIX MUST BE SQUARE", .{});

            for (0..self.matrix.len) |i| for (0..self.matrix.len) |j| {
                rpn_array[i * self.matrix.len + j] = try shuntingYard(T, self.matrix[i][j], self.variables, allocator);
            };

            for (0..self.variables.len) |k| try map.put(self.variables[k], undefined);

            custom_potential_data = .{
                .rpn_array = rpn_array,
                .map = map,
                .allocator = allocator
            };

            return custom_potential_data;
        }
    };
}
