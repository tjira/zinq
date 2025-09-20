//! File with potential that reads values from a file.

const std = @import("std");

const device_read = @import("device_read.zig");
const error_handling = @import("error_handling.zig");
const linear_interpolation = @import("linear_interpolation.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const lerp = linear_interpolation.lerp;
const readRealMatrix = device_read.readRealMatrix;
const throw = error_handling.throw;

/// Struct holding parameters for the file potential.
pub fn FilePotential(comptime T: type) type {
    return struct {
        ndim: u32,
        nstate: u32,
        source : union(enum) {
            path: []const u8,
            data: []T,
        },

        /// Diabatic potential matrix evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), _: T) !void {
            const potential_matrix = try self.getMatrix();

            for (0..U.rows) |i| for (i..U.cols) |j| {
                U.ptr(i, j).* = try lerp(T, potential_matrix, self.ndim + i * U.cols + j, position); U.ptr(j, i).* = U.at(i, j);
            };
        }

        /// Returns the data represented as a matrix. The pointer to the matrix data is stored within this struct.
        pub fn getMatrix(self: @This()) !RealMatrix(T) {
            if (self.source == .path) return throw(RealMatrix(T), "POTENTIAL DATA NOT INITIALIZED, CALL init() FIRST", .{});

            return RealMatrix(T){
                .data = self.source.data,
                .rows = self.source.data.len / (self.nstate * self.nstate + self.ndim),
                .cols = self.nstate * self.nstate + self.ndim,
                .allocator = null
            };
        }

        /// Read the potential data from file, if the file path is specified. Otherwise, just return the matrix from the data pointer.
        pub fn init(self: *@This(), allocator: std.mem.Allocator) !RealMatrix(T) {
            const U = if (self.source == .path) try readRealMatrix(T, self.source.path, allocator) else return self.getMatrix();

            if (self.source == .path) self.source = .{.data = U.data};

            return U;
        }
    };
}
