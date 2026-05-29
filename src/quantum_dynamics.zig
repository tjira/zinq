const std = @import("std");

const matrix = @import("matrix.zig");
const vector = @import("vector.zig");
const writer = @import("writer.zig");

const Matrix = matrix.Matrix;
const Vector = vector.Vector;

const printf = writer.printf;

const InitialConditions = struct {
    position: []const f64,
    momentum: []const f64,
    gamma: []const f64,
    state: u32 = 0,
    adiabatic: bool = false,
};

const Write = struct {
    kinetic_energy: ?[]const u8 = null,
    momentum: ?[]const u8 = null,
    population: ?[]const u8 = null,
    position: ?[]const u8 = null,
    potential_energy: ?[]const u8 = null,
    wavefunction: ?[]const u8 = null,
    total_energy: ?[]const u8 = null,
};

pub const Options = struct {
    grid: struct {
        limits: []const []const f64,
        npoint: u32,
    },
    imaginary: ?struct {
        nstate: u32 = 1,
    } = null,

    initial_conditions: InitialConditions,
    mass: f64,
    iterations: u32,
    time_step: f64,
    write: Write = .{},
    adiabatic: bool = false,
    log_interval: u32 = 1,
};

fn Grid(comptime T: type) type {
    return struct {
        r: Matrix(T),
        k: Matrix(T),
        dr: T,
        dk: T,
    };
}

fn Wavefunction(comptime T: type) type {
    return struct {
        data: Matrix(T),
        decay: Vector(T),

        pub fn init(ndim: u32, nstate: u32, npoint: u32, gpa: std.mem.Allocator) !@This() {
            const rows = std.math.pow(u32, npoint, ndim);

            return @This(){
                .data = try Matrix(T).init(rows, nstate, gpa),
                .decay = try Vector(T).init(nstate, gpa),
            };
        }

        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            self.data.deinit(gpa);
            self.decay.deinit(gpa);
        }
    };
}

pub fn run(comptime _: type, io: std.Io, opt: Options, gpa: std.mem.Allocator) !void {
    var wfn = try Wavefunction(f64).init(1, 2, 128, gpa);
    defer wfn.deinit(gpa);

    _ = opt;

    try printf(io, "{any}\n", .{wfn.data.shape});
}
