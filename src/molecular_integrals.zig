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
const print = device_write.print;
const printJson = device_write.printJson;

/// Integral engine target options.
pub fn Options(comptime _: type) type {
    return struct {
        const Write = struct {
            overlap: ?[]const u8 = null,
            kinetic: ?[]const u8 = null,
            nuclear: ?[]const u8 = null,
            coulomb: ?[]const u8 = null,
        };

        basis: []const u8,
        system: []const u8,

        nthread: u32 = 1,

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
            return @This(){
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
pub fn run(comptime T: type, options: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(options);

    const system = try classical_particle.read(T, options.system, 0, allocator); defer system.deinit();
    var basis = try BasisSet(T).init(system, options.basis, allocator); defer basis.deinit();

    var output = Output(T).init(allocator);

    if (enable_printing) try print("\nNUMBER OF BASIS FUNCTIONS: {d}\n", .{basis.nbf()});

    if (enable_printing and (options.write.overlap != null or options.write.kinetic != null or options.write.nuclear != null or options.write.coulomb != null)) {
        try print("\n", .{});
    }

    if (options.write.overlap) |path| {

        var timer = try std.time.Timer.start();

        if (enable_printing) try print("OVERLAP INTEGRALS: ", .{});

        output.S = try overlap(T, basis, options.nthread, allocator);

        if (enable_printing) try print("{D}\n", .{timer.read()});

        try exportRealMatrix(T, path, output.S.?);
    }

    if (options.write.kinetic) |path| {

        var timer = try std.time.Timer.start();

        if (enable_printing) try print("KINETIC INTEGRALS: ", .{});

        output.K = try kinetic(T, basis, options.nthread, allocator);

        if (enable_printing) try print("{D}\n", .{timer.read()});

        try exportRealMatrix(T, path, output.K.?);
    }

    if (options.write.nuclear) |path| {

        var timer = try std.time.Timer.start();

        if (enable_printing) try print("NUCLEAR INTEGRALS: ", .{});

        output.V = try nuclear(T, system, basis, options.nthread, allocator);

        if (enable_printing) try print("{D}\n", .{timer.read()});

        try exportRealMatrix(T, path, output.V.?);
    }

    if (options.write.coulomb) |path| {

        var timer = try std.time.Timer.start();

        if (enable_printing) try print("COULOMB INTEGRALS: ", .{});

        output.J = try coulomb(T, basis, options.nthread, allocator);

        if (enable_printing) try print("{D}\n", .{timer.read()});

        try exportRealTensorFour(T, path, output.J.?);
    }

    return output;
}

/// Compute the coulomb tensor.
pub fn coulomb(comptime T: type, basis: BasisSet(T), nthread: usize, allocator: std.mem.Allocator) !RealTensor4(T) {
    var J = try RealTensor4(T).init(.{basis.nbf(), basis.nbf(), basis.nbf(), basis.nbf()}, allocator);

    var pool: std.Thread.Pool = undefined; try pool.init(.{.n_jobs = nthread, .allocator = allocator}); defer pool.deinit();

    for (0..J.shape[0]) |i| for (i..J.shape[1]) |j| for (i..J.shape[2]) |k| for ((if (i == k) j else k)..J.shape[3]) |l| {
        try pool.spawn(coulombAssign, .{T, &J, i, j, k, l, basis});
    };

    return J;
}

/// Assigns values to the coulomb tensor based on indices.
pub fn coulombAssign(comptime T: type, J: *RealTensor4(T), i: usize, j: usize, k: usize, l: usize, basis: BasisSet(T)) void {
    const I = basis.contracted_gaussians[i].coulomb(basis.contracted_gaussians[j], basis.contracted_gaussians[k], basis.contracted_gaussians[l]);

    J.ptr(i, j, k, l).* = I;
    J.ptr(i, j, l, k).* = I;
    J.ptr(j, i, k, l).* = I;
    J.ptr(j, i, l, k).* = I;

    J.ptr(k, l, i, j).* = I;
    J.ptr(k, l, j, i).* = I;
    J.ptr(l, k, i, j).* = I;
    J.ptr(l, k, j, i).* = I;
}

/// Compute the kinetic matrix.
pub fn kinetic(comptime T: type, basis: BasisSet(T), nthread: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    var K = try RealMatrix(T).init(basis.nbf(), basis.nbf(), allocator);

    var pool: std.Thread.Pool = undefined; try pool.init(.{.n_jobs = nthread, .allocator = allocator}); defer pool.deinit();

    for (0..K.rows) |i| for (i..K.cols) |j| {
        try pool.spawn(kineticAssign, .{T, &K, i, j, basis});
    };

    return K;
}

/// Assigns values to the kinetic matrix based on indices.
pub fn kineticAssign(comptime T: type, K: *RealMatrix(T), i: usize, j: usize, basis: BasisSet(T)) void {
    K.ptr(i, j).* = basis.contracted_gaussians[i].kinetic(basis.contracted_gaussians[j]); K.ptr(j, i).* = K.at(i, j);
}

/// Compute the nuclear matrix.
pub fn nuclear(comptime T: type, system: ClassicalParticle(T), basis: BasisSet(T), nthread: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    var V = try RealMatrix(T).init(basis.nbf(), basis.nbf(), allocator);

    var pool: std.Thread.Pool = undefined; try pool.init(.{.n_jobs = nthread, .allocator = allocator}); defer pool.deinit();

    for (0..V.rows) |i| for (i..V.cols) |j| {
        try pool.spawn(nuclearAssign, .{T, &V, i, j, system, basis});
    };

    return V;
}

/// Assigns values to the nuclear matrix based on indices.
pub fn nuclearAssign(comptime T: type, V: *RealMatrix(T), i: usize, j: usize, system: ClassicalParticle(T), basis: BasisSet(T)) void {
    V.ptr(i, j).* = basis.contracted_gaussians[i].nuclear(basis.contracted_gaussians[j], system); V.ptr(j, i).* = V.at(i, j);
}

/// Compute the overlap matrix.
pub fn overlap(comptime T: type, basis: BasisSet(T), nthread: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    var S = try RealMatrix(T).init(basis.nbf(), basis.nbf(), allocator);

    var pool: std.Thread.Pool = undefined; try pool.init(.{.n_jobs = nthread, .allocator = allocator}); defer pool.deinit();

    for (0..S.rows) |i| for (i..S.cols) |j| {
        try pool.spawn(overlapAssign, .{T, &S, i, j, basis});
    };

    return S;
}

/// Assigns values to the overlap matrix based on indices.
pub fn overlapAssign(comptime T: type, S: *RealMatrix(T), i: usize, j: usize, basis: BasisSet(T)) void {
    S.ptr(i, j).* = basis.contracted_gaussians[i].overlap(basis.contracted_gaussians[j]); S.ptr(j, i).* = S.at(i, j);
}
