//! General system struct that holds all information about the classical like coordinates, velocities, etc.

const std = @import("std");

const bias_potential = @import("bias_potential.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const string_manipulation = @import("string_manipulation.zig");

const BiasPotential = bias_potential.BiasPotential;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const AN2SM = global_variables.AN2SM;
const throw = error_handling.throw;
const uncr = string_manipulation.uncr;

const WRITE_BUFFER_SIZE = global_variables.WRITE_BUFFER_SIZE;
const A2AU = global_variables.A2AU;
const AN2M = global_variables.AN2M;
const SM2AN = global_variables.SM2AN;
const AU2K = global_variables.AU2K;
const U2AU = global_variables.U2AU;

/// Classical particle struct that holds all information about the classical-like coordinates, velocities and masses.
pub fn ClassicalParticle(comptime T: type) type {
    return struct {
        atoms: ?[]usize,
        acceleration: RealVector(T),
        charge: i32,
        masses: RealVector(T),
        ndim: usize,
        ndof: usize,
        position: RealVector(T),
        velocity: RealVector(T),

        /// Initialize the system using the number of dimensions, an array of masses, and an allocator. The positions and velocities are undefined after initialization.
        pub fn init(ndim: usize, masses: []const T, allocator: std.mem.Allocator) !@This() {
            if (ndim != masses.len) return throw(@This(), "THE LENGTH OF MASSES VECTOR MUST BE THE SAME AS NUMBER OF DIMENSIONS", .{});

            var particle = @This(){
                .atoms = null,
                .acceleration = try RealVector(T).init(ndim, allocator),
                .charge = 0,
                .masses = try RealVector(T).init(ndim, allocator),
                .ndim = ndim,
                .ndof = ndim,
                .position = try RealVector(T).init(ndim, allocator),
                .velocity = try RealVector(T).init(ndim, allocator)
            };

            for (masses, 0..) |mi, i| particle.masses.ptr(i).* = mi;

            return particle;
        }

        /// Initialize the system using the number of dimensions, an array of masses, and an allocator. The positions and velocities are initialized to zero.
        pub fn initZero(ndim: usize, masses: []const T, allocator: std.mem.Allocator) !@This() {
            var particle = try @This().init(ndim, masses, allocator);

            particle.position.zero();
            particle.velocity.zero();
            particle.acceleration.zero();

            return particle;
        }

        /// Free the memory allocated for the system.
        pub fn deinit(self: @This(), allocator: std.mem.Allocator) void {
            if (self.atoms) |atoms| allocator.free(atoms);

            self.position.deinit(allocator);
            self.velocity.deinit(allocator);
            self.acceleration.deinit(allocator);
            self.masses.deinit(allocator);
        }

        /// Calculate the acceleration of the classical particle using the forces from the electronic potential.
        pub fn calculateAcceleration(self: *@This(), pot: ElectronicPotential(T), pot_matrix: *RealMatrix(T), time: T, current_state: usize, fdiff_step: T, bias: ?BiasPotential(T), dir: std.fs.Dir) !void {
            for (0..self.ndim) |i| {

                const force = try pot.forceAdiabatic(pot_matrix, self.position, time, current_state, i, fdiff_step, bias, dir);

                self.acceleration.ptr(i).* = force / self.masses.at(i);
            }
        }

        /// Clone the classical particle.
        pub fn clone(self: @This(), allocator: std.mem.Allocator) !@This() {
            var atoms: ?[]usize = null;

            if (self.atoms != null) {

                atoms = try allocator.alloc(usize, self.atoms.?.len);

                for (0..self.atoms.?.len) |i| atoms.?[i] = self.atoms.?[i];
            }

            return @This(){
                .atoms = atoms,
                .acceleration = try self.acceleration.clone(allocator),
                .charge = self.charge,
                .masses = try self.masses.clone(allocator),
                .ndim = self.ndim,
                .ndof = self.ndof,
                .position = try self.position.clone(allocator),
                .velocity = try self.velocity.clone(allocator)
            };
        }

        /// Calculate the kinetic energy of the classical particle.
        pub fn kineticEnergy(self: @This()) T {
            var kinetic_energy: T = 0;

            for (0..self.ndim) |i| {
                kinetic_energy += 0.5 * self.masses.at(i) * self.velocity.at(i) * self.velocity.at(i);
            }

            return kinetic_energy;
        }

        /// Calculate the kinetic temperature of the system using the equipartition theorem.
        pub fn kineticTemperature(self: @This()) T {
            return 2 * self.kineticEnergy() / @as(T, @floatFromInt(self.ndof)) * AU2K;
        }

        /// Get the number of occupied spatial orbitals.
        pub fn noccSpatial(self: @This()) !usize {
            const spin = try self.noccSpin();

            if (spin % 2 != 0) {
                return throw(usize, "THE NUMBER OF OCCUPIED SPIN ORBITALS MUST BE EVEN TO GET THE NUMBER OF OCCUPIED SPATIAL ORBITALS", .{});
            }

            return spin / 2;
        }

        /// Get the number of occupied spatial orbitals.
        pub fn noccSpin(self: @This()) !usize {
            if (self.atoms == null) return 0;

            var sum: usize = 0;

            for (0..self.atoms.?.len) |i| sum += self.atoms.?[i];

            if (self.charge > @as(i32, @intCast(sum))) {
                return throw(usize, "THE SYSTEM CHARGE CAN'T BE LARGER THAN THE NUMBER OF ELECTRONS", .{});
            }

            return @intCast(@as(i32, @intCast(sum)) - self.charge);
        }

        /// Get the nuclear repulsion energy of the classical particle.
        pub fn nuclearRepulsionEnergy(self: @This()) T {
            if (self.atoms == null) return 0;

            var energy: T = 0;

            for (0..self.atoms.?.len) |i| for (i + 1..self.atoms.?.len) |j| {

                const dx = self.position.at(3 * i + 0) - self.position.at(3 * j + 0);
                const dy = self.position.at(3 * i + 1) - self.position.at(3 * j + 1);
                const dz = self.position.at(3 * i + 2) - self.position.at(3 * j + 2);

                const r = std.math.sqrt(dx * dx + dy * dy + dz * dz);

                energy += @as(T, @floatFromInt(self.atoms.?[i] * self.atoms.?[j])) / r;
            };

            return energy;
        }

        /// First half of the velocity verlet algorithm, which updates the positions and velocities using the current accelerations. The second half of the algorithm should be called after updating the accelerations.
        pub fn propagateVelocityVerletFirstHalf(self: *@This(), time_step: T) void {
            for (0..self.ndim) |i| {
                self.velocity.ptr(i).* += 0.5 * self.acceleration.at(i) * time_step;
                self.position.ptr(i).* += self.velocity.at(i) * time_step;
            }
        }

        /// Second half of the velocity verlet algorithm, which updates the velocities using the new accelerations. The first half of the algorithm should be called before updating the accelerations.
        pub fn propagateVelocityVerletSecondHalf(self: *@This(), time_step: T) !void {
            for (0..self.ndim) |i| {
                self.velocity.ptr(i).* += 0.5 * self.acceleration.at(i) * time_step;
            }
        }

        /// Sets the position of the particle from normal distribution with given mean and standard deviation using the provided random number generator state.
        pub fn setPositionRandn(self: *@This(), mean: []const T, stdev: []const T, random: *std.Random) !void {
            if (mean.len != self.ndim or stdev.len != self.ndim) return throw(void, "POSITION MEAN AND STDEV MUST HAVE THE SAME DIMENSION AS THE SYSTEM", .{});

            for (0..self.ndim) |i| {
                self.position.ptr(i).* = mean[i] + stdev[i] * random.floatNorm(T);
            }
        }

        /// Sets the momentum of the particle from normal distribution with given mean and standard deviation using the provided random number generator state.
        pub fn setMomentumRandn(self: *@This(), mean: []const T, stdev: []const T, random: *std.Random) !void {
            if (mean.len != self.ndim or stdev.len != self.ndim) return throw(void, "MOMENTUM MEAN AND STDEV MUST HAVE THE SAME DIMENSION AS THE SYSTEM", .{});

            for (0..self.ndim) |i| {
                self.velocity.ptr(i).* = (mean[i] + stdev[i] * random.floatNorm(T)) / self.masses.at(i);
            }
        }

        /// Write the coordinates of the system to an .xyz file. The coordinates are converted from atomic units to Angstroms before writing.
        pub fn writeCoordinatesToXYZ(self: @This(), filename: []const u8, dir: std.fs.Dir) !void {
            if (self.atoms == null) return throw(void, "CAN'T WRITE COORDINATES TO FILE IF ATOMS VECTOR IS NULL", .{});

            const file = try dir.createFile(filename, .{}); defer file.close();

            var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

            var writer = file.writer(&buffer); var writer_interface = &writer.interface;

            try writer_interface.print("{d}\n\n", .{self.position.len / 3});

            for (0..self.position.len / 3) |i| {

                const x = self.position.at(i * 3 + 0) / A2AU;
                const y = self.position.at(i * 3 + 1) / A2AU;
                const z = self.position.at(i * 3 + 2) / A2AU;
                
                try writer_interface.print("{s:2} {d:20.14} {d:20.14} {d:20.14}\n", .{try AN2SM(self.atoms.?[i]), x, y, z});
            }

            try writer_interface.flush();
        }
    };
}

/// Read the system from a .xyz file.
pub fn read(comptime T: type, path: []const u8, charge: i32, geometry: usize, allocator: std.mem.Allocator) !ClassicalParticle(T) {
    const file = try std.fs.cwd().openFile(path, .{}); defer file.close();

    var buffer: [1024]u8 = undefined; var reader = file.reader(&buffer); var reader_interface = &reader.interface;

    const natom = try std.fmt.parseInt(u32, uncr(try reader_interface.takeDelimiterExclusive('\n')), 10);

    var atoms = try allocator.alloc(usize, natom);

    reader_interface.toss(1); _ = try reader_interface.discardDelimiterInclusive('\n');

    for (0..(natom + 2) * geometry) |_| _ = try reader_interface.discardDelimiterInclusive('\n');

    var position = try RealVector(T).init(3 * natom, allocator); defer position.deinit(allocator);

    for (0..natom) |i| {

        var it = std.mem.tokenizeAny(u8, try reader_interface.takeDelimiterExclusive('\n'), " "); 

        atoms[i] = SM2AN.get(it.next().?).?;

        for (0..3) |j| {
            position.ptr(3 * i + j).* = try std.fmt.parseFloat(T, uncr(it.next().?)) * A2AU;
        }

        reader_interface.toss(1);
    }

    var masses = try allocator.alloc(T, 3 * natom); defer allocator.free(masses);

    for (0..3 * natom) |i| masses[i] = AN2M[atoms[i / 3]] * U2AU;

    var system = try ClassicalParticle(T).initZero(3 * natom, masses, allocator);

    system.atoms = atoms; system.charge = charge; system.ndof -= 6; try position.copyTo(&system.position);

    return system;
}
