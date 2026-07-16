//! Command-line entry point for running molecular dynamics and quantum chemistry simulations.

const builtin = @import("builtin");
const config = @import("config");
const std = @import("std");

const libint = @cImport(@cInclude("libint2/config.h"));
const fftw = @cImport(@cInclude("fftw3.h"));
const libxc = @cImport(@cInclude("xc.h"));
const openblas = @cImport(@cInclude("openblas_config.h"));

const Allocator = std.mem.Allocator;

pub const cimport = @import("cimport.zig");
pub const classical_dynamics = @import("classical_dynamics.zig");
pub const configuration_interaction = @import("configuration_interaction.zig");
pub const constant = @import("constant.zig");
pub const cphf = @import("cphf.zig");
pub const density_functional_theory = @import("density_functional_theory.zig");
pub const dual = @import("dual.zig");
pub const ehrenfest = @import("ehrenfest.zig");
pub const exprtk = @import("exprtk.zig");
pub const fourier_transform = @import("fourier_transform.zig");
pub const frequency_analysis = @import("frequency_analysis.zig");
pub const hartree_fock = @import("hartree_fock.zig");
pub const integral_transform = @import("integral_transform.zig");
pub const integrator = @import("integrator.zig");
pub const lebedev_quadrature_nodes = @import("lebedev_quadrature_nodes.zig");
pub const linear_algebra = @import("linear_algebra.zig");
pub const molecular_grid = @import("molecular_grid.zig");
pub const molecular_integrals = @import("molecular_integrals.zig");
pub const molecular_optimization = @import("molecular_optimization.zig");
pub const molecular_system = @import("molecular_system.zig");
pub const moller_plesset = @import("moller_plesset.zig");
pub const nuclear_derivative = @import("nuclear_derivative.zig");
pub const parser = @import("parser.zig");
pub const population_analysis = @import("population_analysis.zig");
pub const potential = @import("potential.zig");
pub const potential_plot = @import("potential_plot.zig");
pub const quantum_dynamics = @import("quantum_dynamics.zig");
pub const read_write = @import("read_write.zig");
pub const spectral_analysis = @import("spectral_analysis.zig");
pub const surface_hopping = @import("surface_hopping.zig");
pub const tensor = @import("tensor.zig");
pub const value = @import("value.zig");
pub const wavepacket = @import("wavepacket.zig");
pub const xc_functional = @import("xc_functional.zig");

const printf = read_write.printf;

/// Simulation runner handlers for classical dynamics, Hartree-Fock, and quantum dynamics.
const Handlers = struct {
    pub const classical_dynamics = @import("classical_dynamics.zig");
    pub const configuration_interaction = @import("configuration_interaction.zig");
    pub const hartree_fock = @import("hartree_fock.zig");
    pub const molecular_integrals = @import("molecular_integrals.zig");
    pub const moller_plesset = @import("moller_plesset.zig");
    pub const potential_plot = @import("potential_plot.zig");
    pub const quantum_dynamics = @import("quantum_dynamics.zig");
};

/// Parsed configuration options representing target molecular dynamics or electronic structure jobs.
const Options = struct {
    zinq: []union(enum) {
        classical_dynamics: classical_dynamics.Options,
        configuration_interaction: configuration_interaction.Options,
        hartree_fock: hartree_fock.Options,
        molecular_integrals: molecular_integrals.Options,
        moller_plesset: moller_plesset.Options,
        potential_plot: potential_plot.Options,
        quantum_dynamics: quantum_dynamics.Options,
    },
};

/// Main entry point printing library versions and executing molecular simulation targets.
pub fn main(init: std.process.Init) !void {
    var timer = std.Io.Timestamp.now(init.io, .real);

    const v_major = builtin.zig_version.major;
    const v_minor = builtin.zig_version.minor;
    const v_patch = builtin.zig_version.patch;

    try printf(init.io, "ZIG: v{d}.{d}.{d}, ZINQ: {s}\n\n", .{ v_major, v_minor, v_patch, config.version });

    const openblas_v = std.mem.trim(u8, openblas.OPENBLAS_VERSION, "OpenBLAS ");

    try printf(init.io, "OPENBLAS: v{s}, ", .{openblas_v});

    const fftw_v = std.mem.trim(u8, std.mem.span(fftw.fftw_version), "fftw-");

    try printf(init.io, "FFTW: v{s}, ", .{fftw_v});

    try printf(init.io, "LIBINT: v{s}, LIBXC: v{s}\n", .{ libint.LIBINT_VERSION, libxc.xc_version_string() });

    const args = try init.minimal.args.toSlice(init.arena.allocator());

    var parsobj = try parser.Parser.init(args, init.arena.allocator());

    try parsobj.dispatch(init.io, init.gpa, init.arena.allocator());

    if (parsobj.action == .help) return;

    try printf(init.io, "\nTOTAL EXECUTION TIME: {f}\n", .{timer.untilNow(init.io, .real)});
}

/// Reads and parses the JSON configuration file containing molecular simulation parameters.
fn parse(comptime T: type, io: std.Io, fname: []const u8, arena: Allocator) !?std.json.Parsed(T) {
    const fcontent = std.Io.Dir.cwd().readFileAlloc(io, fname, arena, .unlimited) catch |err| {
        if (err == error.FileNotFound) {
            try printf(io, "\nINPUT FILE '{s}' NOT FOUND\n", .{fname});

            return null;
        }

        return err;
    };

    return try std.json.parseFromSlice(T, arena, fcontent, .{});
}

/// Runs the specified electronic structure or molecular dynamics jobs.
pub fn run(comptime T: type, io: std.Io, fname: []const u8, gpa: Allocator, arena: Allocator) !void {
    const parsed = try parse(Options, io, fname, arena) orelse return;

    for (0..parsed.value.zinq.len) |i| {
        try printf(io, "\nRUNNING TARGET: {s}/#{d}\n", .{ fname, i + 1 });

        switch (parsed.value.zinq[i]) {
            inline else => |field, tag| {
                var result = try @field(Handlers, @tagName(tag)).run(T, io, field, true, gpa);
                defer result.deinit(gpa);
            },
        }
    }
}
