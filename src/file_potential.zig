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

var file_potential_data: ?FilePotentialData(f64) = null;

/// Struct holding the data for the file-based potential.
pub fn FilePotentialData(comptime T: type) type {
    return struct {
        data: RealMatrix(T),
        allocator: std.mem.Allocator,

        /// Free the resources.
        pub fn deinit(self: *@This()) void {
            self.data.deinit();
        }
    };
}

/// Struct holding parameters for the file potential.
pub fn FilePotential(comptime T: type) type {
    return struct {
        ndim: u32,
        nstate: u32,
        path: []const u8,

        /// Diabatic potential matrix evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), time: T) !void {
            for (0..U.rows) |i| for (i..U.cols) |j| {
                U.ptr(i, j).* = try evaluateDiabaticElement(self, i, j, position, time); U.ptr(j, i).* = U.at(i, j);
            };
        }

        /// Diabatic potential matrix element evaluator.
        pub fn evaluateDiabaticElement(self: @This(), i: usize, j: usize, position: RealVector(T), _: T) !T {
            if (file_potential_data == null) return throw(T, "POTENTIAL DATA NOT INITIALIZED, CALL init() FIRST", .{});
            if (i >= self.nstate or j >= self.nstate) return throw(T, "INVALID INDEX WHEN EVALUATING DIABATIC MATRIX ELEMENT", .{});

            return try lerp(T, file_potential_data.?.data, self.ndim + i * self.nstate + j, position);
        }

        /// Read the potential data from file, if the file path is specified. Otherwise, just return the matrix from the data pointer.
        pub fn init(self: @This(), allocator: std.mem.Allocator) !?FilePotentialData(T) {
            const U = try readRealMatrix(T, self.path, allocator);

            if (self.ndim + self.nstate * self.nstate != U.cols) return throw(?FilePotentialData(T), "POTENTIAL DATA DIMENSIONS DO NOT MATCH THE SPECIFIED NDIM AND NSTATE", .{});

            file_potential_data = .{
                .data = U,
                .allocator = allocator
            };

            return file_potential_data;
        }
    };
}
