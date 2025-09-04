//! Integral engine.

const std = @import("std");

const basis_set = @import("basis_set.zig");
const classical_particle = @import("classical_particle.zig");
const contracted_gaussian = @import("contracted_gaussian.zig");
const device_write = @import("device_write.zig");
const real_matrix = @import("real_matrix.zig");
const real_tensor_four = @import("real_tensor_four.zig");

const BasisSet = basis_set.BasisSet;
const ClassicalParticle = classical_particle.ClassicalParticle;
const ContractedGaussian = contracted_gaussian.ContractedGaussian;
const RealMatrix = real_matrix.RealMatrix;
const RealTensor4 = real_tensor_four.RealTensor4;
const exportRealMatrix = device_write.exportRealMatrix;
const exportRealTensorFour = device_write.exportRealTensorFour;

/// Integral engine target oprions.
pub fn Options(comptime _: type) type {
    return struct {
        const Write = struct {
            overlap: ?[]const u8 = null,
            kinetic: ?[]const u8 = null,
            nuclear: ?[]const u8 = null,
            coulomb: ?[]const u8 = null,
        };

        system: []const u8,
        basis: []const u8,

        write: Write = .{}
    };
}

/// Integral engine target output.
pub fn Output(comptime T: type) type {
    return struct {
        S: ?RealMatrix(T),
        K: ?RealMatrix(T),
        V: ?RealMatrix(T),
        J: ?RealTensor4(T),

        allocator: std.mem.Allocator,

        /// Initialize the output struct.
        pub fn init(allocator: std.mem.Allocator) @This() {
            return Output(T){
                .S = null,
                .K = null,
                .V = null,
                .J = null,
                .allocator = allocator,
            };
        }

        /// Deinitialize the output struct.
        pub fn deinit(self: @This()) void {
            if (self.S) |S| S.deinit();
            if (self.K) |K| K.deinit();
            if (self.V) |V| V.deinit();
            if (self.J) |J| J.deinit();
        }
    };
}

/// Run the integral engine target.
pub fn run(comptime T: type, options: Options(T), _: bool, allocator: std.mem.Allocator) !Output(T) {
    const system = try classical_particle.read(T, options.system, 0, allocator); defer system.deinit();
    var basis = try BasisSet(T).init(system, options.basis, allocator); defer basis.deinit();

    var output = Output(T).init(allocator);

    if (options.write.overlap) |path| {

        output.S = try overlap(T, basis, allocator);

        try exportRealMatrix(T, path, output.S.?);
    }

    if (options.write.kinetic) |path| {

        output.K = try kinetic(T, basis, allocator);

        try exportRealMatrix(T, path, output.K.?);
    }

    if (options.write.nuclear) |path| {

        output.V = try nuclear(T, system, basis, allocator);

        try exportRealMatrix(T, path, output.V.?);
    }

    if (options.write.coulomb) |path| {

        output.J = try coulomb(T, basis, allocator);

        try exportRealTensorFour(T, path, output.J.?);
    }

    return output;
}

/// Compute the coulomb tensor.
pub fn coulomb(comptime T: type, basis: BasisSet(T), allocator: std.mem.Allocator) !RealTensor4(T) {
    var J = try RealTensor4(T).init(.{basis.contracted_gaussians.len, basis.contracted_gaussians.len, basis.contracted_gaussians.len, basis.contracted_gaussians.len}, allocator);

    for (0..J.shape[0]) |i| for (i..J.shape[1]) |j| for (i..J.shape[2]) |k| for ((if (i == k) j else k)..J.shape[3]) |l| {

        const I = basis.contracted_gaussians[i].coulomb(basis.contracted_gaussians[j], basis.contracted_gaussians[k], basis.contracted_gaussians[l]);

        J.ptr(i, j, k, l).* = I;
        J.ptr(i, j, l, k).* = I;
        J.ptr(j, i, k, l).* = I;
        J.ptr(j, i, l, k).* = I;

        J.ptr(k, l, i, j).* = I;
        J.ptr(k, l, j, i).* = I;
        J.ptr(l, k, i, j).* = I;
        J.ptr(l, k, j, i).* = I;
    };

    return J;
}

/// Compute the kinetic matrix.
pub fn kinetic(comptime T: type, basis: BasisSet(T), allocator: std.mem.Allocator) !RealMatrix(T) {
    var K = try RealMatrix(T).init(basis.contracted_gaussians.len, basis.contracted_gaussians.len, allocator);

    for (0..K.rows) |i| for (i..K.cols) |j| {
        K.ptr(i, j).* = basis.contracted_gaussians[i].kinetic(basis.contracted_gaussians[j]); K.ptr(j, i).* = K.at(i, j);
    };

    return K;
}

/// Compute the kinetic matrix.
pub fn nuclear(comptime T: type, system: ClassicalParticle(T), basis: BasisSet(T), allocator: std.mem.Allocator) !RealMatrix(T) {
    var V = try RealMatrix(T).init(basis.contracted_gaussians.len, basis.contracted_gaussians.len, allocator);

    for (0.. V.rows) |i| for (i..V.cols) |j| {
        V.ptr(i, j).* = basis.contracted_gaussians[i].nuclear(basis.contracted_gaussians[j], system); V.ptr(j, i).* = V.at(i, j);
    };

    return V;
}

/// Compute the overlap matrix.
pub fn overlap(comptime T: type, basis: BasisSet(T), allocator: std.mem.Allocator) !RealMatrix(T) {
    var S = try RealMatrix(T).init(basis.contracted_gaussians.len, basis.contracted_gaussians.len, allocator);

    for (0..S.rows) |i| for (i..S.cols) |j| {
        S.ptr(i, j).* = basis.contracted_gaussians[i].overlap(basis.contracted_gaussians[j]); S.ptr(j, i).* = S.at(i, j);
    };

    return S;
}
