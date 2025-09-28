//! Integral engine.

const std = @import("std");

const basis_set = @import("basis_set.zig");
const classical_particle = @import("classical_particle.zig");
const contracted_gaussian = @import("contracted_gaussian.zig");
const device_write = @import("device_write.zig");
const parallel_algorithm = @import("parallel_algorithm.zig");
const real_matrix = @import("real_matrix.zig");
const real_tensor_four = @import("real_tensor_four.zig");

const BasisSet = basis_set.BasisSet;
const ClassicalParticle = classical_particle.ClassicalParticle;
const ContractedGaussian = contracted_gaussian.ContractedGaussian;
const RealMatrix = real_matrix.RealMatrix;
const RealTensor4 = real_tensor_four.RealTensor4;

const exportRealMatrix = device_write.exportRealMatrix;
const exportRealTensorFour = device_write.exportRealTensorFour;
const parallelFor = parallel_algorithm.parallelFor;
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

    const func = struct {
        fn call(i: usize, args: anytype) void {

            for (i..args.J.shape[1]) |j| for (i..args.J.shape[2]) |k| for ((if (i == k) j else k)..args.J.shape[3]) |l| {

                const I = args.basis.contracted_gaussians[i].coulomb(args.basis.contracted_gaussians[j], args.basis.contracted_gaussians[k], args.basis.contracted_gaussians[l]);

                args.J.ptr(i, j, k, l).* = I;
                args.J.ptr(i, j, l, k).* = I;
                args.J.ptr(j, i, k, l).* = I;
                args.J.ptr(j, i, l, k).* = I;

                args.J.ptr(k, l, i, j).* = I;
                args.J.ptr(k, l, j, i).* = I;
                args.J.ptr(l, k, i, j).* = I;
                args.J.ptr(l, k, j, i).* = I;
            };
        }
    }.call;

    try parallelFor(func, .{.J = &J, .basis = basis}, 0, J.shape[0], nthread, allocator);

    return J;
}

/// Compute the kinetic matrix.
pub fn kinetic(comptime T: type, basis: BasisSet(T), nthread: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    var K = try RealMatrix(T).init(basis.nbf(), basis.nbf(), allocator);

    const func = struct {
        fn call(i: usize, args: anytype) void {
            for (i..args.K.cols) |j| {
                args.K.ptr(i, j).* = args.basis.contracted_gaussians[i].kinetic(args.basis.contracted_gaussians[j]); args.K.ptr(j, i).* = args.K.at(i, j);
            }
        }
    }.call;

    try parallelFor(func, .{.K = &K, .basis = basis}, 0, K.rows, nthread, allocator);

    return K;
}

/// Compute the kinetic matrix.
pub fn nuclear(comptime T: type, system: ClassicalParticle(T), basis: BasisSet(T), nthread: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    var V = try RealMatrix(T).init(basis.nbf(), basis.nbf(), allocator);

    const func = struct {
        fn call(i: usize, args: anytype) void {
            for (i..args.V.cols) |j| {
                args.V.ptr(i, j).* = args.basis.contracted_gaussians[i].nuclear(args.basis.contracted_gaussians[j], args.system); args.V.ptr(j, i).* = args.V.at(i, j);
            }
        }
    }.call;

    try parallelFor(func, .{.V = &V, .basis = basis, .system = system}, 0, V.rows, nthread, allocator);

    return V;
}

/// Compute the overlap matrix.
pub fn overlap(comptime T: type, basis: BasisSet(T), nthread: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    var S = try RealMatrix(T).init(basis.nbf(), basis.nbf(), allocator);

    const func = struct {
        fn call(i: usize, args: anytype) void {
            for (i..args.S.cols) |j| {
                args.S.ptr(i, j).* = args.basis.contracted_gaussians[i].overlap(args.basis.contracted_gaussians[j]); args.S.ptr(j, i).* = args.S.at(i, j);
            }
        }
    }.call;

    try parallelFor(func, .{.S = &S, .basis = basis}, 0, S.rows, nthread, allocator);

    return S;
}
