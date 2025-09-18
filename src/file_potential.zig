//! File with potential that reads values from a file.

const std = @import("std");

const device_read = @import("device_read.zig");
const device_write = @import("device_write.zig");
const linear_interpolation = @import("linear_interpolation.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const lerp = linear_interpolation.lerp;
const printRealMatrix = device_write.printRealMatrix;
const readRealMatrix = device_read.readRealMatrix;

/// Struct holding parameters for the file potential.
pub fn FilePotential(comptime T: type) type {
    return struct {
        ndim: u32, nstate: u32, path: ?[]const u8, data: ?[]T = null,

        /// Diabatic potential matrix evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), time: T) !void {
            _ = time;

            const potential_matrix = try self.getMatrix();

            for (0..U.rows) |i| for (i..U.cols) |j| {
                U.ptr(i, j).* = try lerp(T, potential_matrix, self.ndim + i * U.cols + j, position); U.ptr(j, i).* = U.at(i, j);
            };
        }

        /// Returns the data represented as a matrix. The pointer to the matrix data is stored within this struct.
        pub fn getMatrix(self: @This()) !RealMatrix(T) {
            if (self.data == null) return error.DataNotRead;

            return RealMatrix(T){
                .data = self.data.?,
                .rows = self.data.?.len / (self.nstate * self.nstate + self.ndim),
                .cols = self.nstate * self.nstate + self.ndim,
                .allocator = null
            };
        }

        /// Read the potential data from file. The function returns the read matrix and the pointer to the matrix data is stored within this struct.
        pub fn read(self: *@This(), allocator: std.mem.Allocator) !RealMatrix(T) {
            const U = try readRealMatrix(T, self.path.?, allocator);

            self.data = U.data;

            return U;
        }
    };
}
