//! File with potential that reads values from a file.

const std = @import("std");

const device_read = @import("device_read.zig");
const linear_interpolation = @import("linear_interpolation.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const lerp = linear_interpolation.lerp;
const readRealMatrix = device_read.readRealMatrix;

/// Struct holding parameters for the file potential.
pub fn FilePotential(comptime T: type) type {
    return struct {
        ndim: u32,
        nstate: u32,
        path: []const u8,

        data: ?struct {
            matrix: RealMatrix(T),
        } = null,

        /// Diabatic potential matrix evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), time: T) !void {
            for (0..U.rows) |i| for (i..U.cols) |j| {
                U.ptr(i, j).* = try evaluateDiabaticElement(self, i, j, position, time);
                U.ptr(j, i).* = U.at(i, j);
            };
        }

        /// Diabatic potential matrix element evaluator.
        pub fn evaluateDiabaticElement(self: @This(), i: usize, j: usize, position: RealVector(T), _: T) !T {
            if (self.data == null) {

                std.log.err("FILE POTENTIAL DATA IS NOT INITIALIZED, MAKE SURE TO CALL init() BEFORE EVALUATING THE FILE POTENTIAL", .{});

                return error.ProgrammingError;
            }

            if (i >= self.nstate or j >= self.nstate) {

                std.log.err("INVALID INDEX ({d}, {d}) WHEN EVALUATING DIABATIC MATRIX ELEMENT, THE POTENTIAL MATRIX IS {d}X{d}", .{i, j, self.nstate, self.nstate});

                return error.ProgrammingError;
            }

            return try lerp(T, self.data.?.matrix, self.ndim + i * self.nstate + j, position);
        }

        /// Read the potential data from file, if the file path is specified. Otherwise, just return the matrix from the data pointer.
        pub fn init(self: *@This(), allocator: std.mem.Allocator) !void {
            const U = try readRealMatrix(T, self.path, allocator);

            if (self.ndim + self.nstate * self.nstate != U.cols) {

                std.log.err("INVALID NUMBER OF COLUMNS IN THE FILE, EXPECTED {d} BUT GOT {d}", .{self.ndim + self.nstate * self.nstate, U.cols});

                return error.InputError;
            }

            self.data = .{
                .matrix = U
            };
        }

        /// Free the resources.
        pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
            if (self.data) |*data| data.matrix.deinit(allocator);
        }

        /// Custom potential JSON parser.
        pub fn jsonParseFromValue(allocator: std.mem.Allocator, source: std.json.Value, options: std.json.ParseOptions) !@This() {
            if (source != .object) return error.UnexpectedToken;

            const ndim_val = source.object.get("ndim") orelse return error.MissingField;
            const nstate_val = source.object.get("nstate") orelse return error.MissingField;
            const path_val = source.object.get("path") orelse return error.MissingField;

            const parsed_ndim = try std.json.innerParseFromValue(u32, allocator, ndim_val, options);
            const parsed_nstate = try std.json.innerParseFromValue(u32, allocator, nstate_val, options);
            const parsed_path = try std.json.innerParseFromValue([]const u8, allocator, path_val, options);

            return .{.ndim = parsed_ndim, .nstate = parsed_nstate, .path = parsed_path};
        }

        /// Custom potential JSON stringifier.
        pub fn jsonStringify(self: *const @This(), jws: anytype) !void {
            try jws.beginObject();

            try jws.objectField("ndim");
            try jws.write(self.ndim);

            try jws.objectField("nstate");
            try jws.write(self.nstate);

            try jws.objectField("path");
            try jws.write(self.path);
            
            try jws.endObject();
        }
    };
}
