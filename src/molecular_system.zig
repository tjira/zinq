//! Represents a molecular system, managing its nuclear geometry, basis set information, and quantum chemical operators.

const std = @import("std");

const libint = @import("cimport.zig").libint;

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

/// Generates a MolecularSystem struct type for the given floating-point coordinate precision type T.
pub fn MolecularSystem(comptime T: type) type {
    // Molecular system containing atom types, 3D coordinates, basis functions, and electron count.
    return struct {
        ptr: *libint.SystemData,

        coors: []T,

        atoms: []i32,
        bf2at: []i32,

        nbf: usize,
        nel: usize,

        /// Initializes the molecular system structure from geometry and basis inputs, calculating the net electron count.
        pub fn init(system: []const u8, basis: []const u8, charge: i32, multiplicity: u32, gpa: Allocator) !@This() {
            const sys_c = try gpa.dupeSentinel(u8, system, 0);
            defer gpa.free(sys_c);

            const bas_c = try gpa.dupeSentinel(u8, basis, 0);
            defer gpa.free(bas_c);

            const ptr = libint.libint_init(sys_c.ptr, bas_c.ptr) orelse return error.InitializationFailed;
            errdefer libint.libint_deinit(ptr);

            const nat = libint.libint_nat(ptr);
            const nbf = libint.libint_nbf(ptr);

            const atoms = try gpa.alloc(i32, 1 * nat);
            errdefer gpa.free(atoms);

            const coors = try gpa.alloc(T, 3 * nat);
            errdefer gpa.free(coors);

            const bf2at = try gpa.alloc(i32, nbf);
            errdefer gpa.free(bf2at);

            libint.libint_atoms(atoms.ptr, ptr);
            libint.libint_coors(coors.ptr, ptr);
            libint.libint_bf2at(bf2at.ptr, ptr);

            var nel: usize = 0;

            for (0..atoms.len) |i| {
                nel += @intCast(atoms[i]);
            }

            if (@as(i32, @intCast(nel)) - charge <= 0) {
                std.log.err("NUMBER OF ELECTRONS MUST BE GREATER THAN 0", .{});

                return error.InvalidInput;
            }

            nel = @intCast(@as(i32, @intCast(nel)) - charge);

            if ((nel + multiplicity - 1) % 2 != 0) {
                std.log.err("CHARGE AND MULTIPLICITY ARE INCOMPATIBLE WITH THE NUMBER OF ELECTRONS", .{});

                return error.InvalidInput;
            }

            if (nel + 1 < multiplicity) {
                std.log.err("MULTIPLICITY IS INCOMPATIBLE WITH THE NUMBER OF ELECTRONS", .{});

                return error.InvalidInput;
            }

            return .{ .ptr = ptr, .nbf = nbf, .atoms = atoms, .coors = coors, .bf2at = bf2at, .nel = nel };
        }

        /// Frees allocated memory for atoms, coordinates, basis functions, and the underlying Libint system data.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            gpa.free(self.atoms);
            gpa.free(self.coors);
            gpa.free(self.bf2at);

            libint.libint_deinit(self.ptr);
        }

        /// Calculates the four-center, two-electron Coulomb repulsion integrals over the basis functions.
        pub fn coulomb(self: @This(), gpa: Allocator) !Tensor(T, 4) {
            const I = try Tensor(T, 4).initZero(.{ self.nbf, self.nbf, self.nbf, self.nbf }, gpa);
            errdefer I.deinit(gpa);

            libint.libint_coulomb(I.data.ptr, self.ptr);

            return I;
        }

        /// Computes the first-order derivatives of the two-electron Coulomb integrals with respect to nuclear coordinates.
        pub fn coulombD1(self: @This(), gpa: Allocator) !Tensor(T, 5) {
            const shape = .{ 3 * libint.libint_nat(self.ptr), self.nbf, self.nbf, self.nbf, self.nbf };

            const I = try Tensor(T, 5).initZero(shape, gpa);
            errdefer I.deinit(gpa);

            libint.libint_coulomb_deriv(I.data.ptr, self.ptr);

            return I;
        }

        /// Computes the first-order derivatives of spin-blocked two-electron Coulomb integrals.
        pub fn coulombD1Spin(self: @This(), gpa: Allocator) !Tensor(T, 5) {
            var J = try self.coulombD1(gpa);
            defer J.deinit(gpa);

            const shape = .{ 3 * libint.libint_nat(self.ptr), 2 * self.nbf, 2 * self.nbf, 2 * self.nbf, 2 * self.nbf };

            var I = try Tensor(T, 5).initZero(shape, gpa);
            errdefer I.deinit(gpa);

            for (0..I.shape[0]) |p| for (0..self.nbf) |i| for (0..self.nbf) |j| for (0..self.nbf) |k| for (0..self.nbf) |l| {
                const val = J.at(.{ p, i, j, k, l });

                I.ptr(.{ p, i, j, k, l }).* = val;

                I.ptr(.{ p, i + self.nbf, j, k + self.nbf, l }).* = val;
                I.ptr(.{ p, i, j + self.nbf, k, l + self.nbf }).* = val;

                I.ptr(.{ p, i + self.nbf, j + self.nbf, k + self.nbf, l + self.nbf }).* = val;
            };

            return I;
        }

        /// Constructs the spin-blocked four-center two-electron Coulomb repulsion integral tensor.
        pub fn coulombSpin(self: @This(), gpa: Allocator) !Tensor(T, 4) {
            var J = try self.coulomb(gpa);
            defer J.deinit(gpa);

            const shape = .{ 2 * self.nbf, 2 * self.nbf, 2 * self.nbf, 2 * self.nbf };

            var I = try Tensor(T, 4).initZero(shape, gpa);
            errdefer I.deinit(gpa);

            for (0..self.nbf) |i| for (0..self.nbf) |j| for (0..self.nbf) |k| for (0..self.nbf) |l| {
                const val = J.at(.{ i, j, k, l });

                I.ptr(.{ i, j, k, l }).* = val;

                I.ptr(.{ i + self.nbf, j, k + self.nbf, l }).* = val;
                I.ptr(.{ i, j + self.nbf, k, l + self.nbf }).* = val;

                I.ptr(.{ i + self.nbf, j + self.nbf, k + self.nbf, l + self.nbf }).* = val;
            };

            return I;
        }

        /// Computes the one-electron kinetic energy matrix elements in the molecular basis.
        pub fn kinetic(self: @This(), gpa: Allocator) !Matrix(T) {
            const I = try Matrix(T).initZero(self.nbf, self.nbf, gpa);
            errdefer I.deinit(gpa);

            libint.libint_kinetic(I.data.ptr, self.ptr);

            return I;
        }

        /// Computes the first-order derivative of the one-electron kinetic energy matrix with respect to nuclear coordinates.
        pub fn kineticD1(self: @This(), gpa: Allocator) !Tensor(T, 3) {
            const shape = .{ 3 * libint.libint_nat(self.ptr), self.nbf, self.nbf };

            const I = try Tensor(T, 3).initZero(shape, gpa);
            errdefer I.deinit(gpa);

            libint.libint_kinetic_deriv(I.data.ptr, self.ptr);

            return I;
        }

        /// Computes the spin-blocked first-order derivative of the kinetic energy matrix.
        pub fn kineticD1Spin(self: @This(), gpa: Allocator) !Tensor(T, 3) {
            var K = try self.kineticD1(gpa);
            defer K.deinit(gpa);

            const shape = .{ 3 * libint.libint_nat(self.ptr), 2 * self.nbf, 2 * self.nbf };

            var I = try Tensor(T, 3).initZero(shape, gpa);
            errdefer I.deinit(gpa);

            for (0..I.shape[0]) |k| for (0..self.nbf) |i| for (0..self.nbf) |j| {
                I.ptr(.{ k, i, j }).* = K.at(.{ k, i, j });

                I.ptr(.{ k, i + self.nbf, j + self.nbf }).* = I.at(.{ k, i, j });
            };

            return I;
        }

        /// Constructs the spin-blocked one-electron kinetic energy matrix.
        pub fn kineticSpin(self: @This(), gpa: Allocator) !Matrix(T) {
            var K = try self.kinetic(gpa);
            defer K.deinit(gpa);

            var I = try Matrix(T).initZero(2 * self.nbf, 2 * self.nbf, gpa);
            errdefer I.deinit(gpa);

            for (0..self.nbf) |i| for (0..self.nbf) |j| {
                I.ptr(i, j).* = K.at(i, j);

                I.ptr(i + self.nbf, j + self.nbf).* = I.at(i, j);
            };

            return I;
        }

        /// Calculates the classical nuclear-nuclear electrostatic repulsion energy of the molecular system.
        pub fn nrep(self: @This()) !T {
            if (self.atoms.len == 0) return 0;

            var energy: T = 0;

            for (0..self.atoms.len) |i| {
                const Zi = @as(T, @floatFromInt(self.atoms[i]));

                const xi = self.coors[3 * i + 0];
                const yi = self.coors[3 * i + 1];
                const zi = self.coors[3 * i + 2];

                for (0..i) |j| {
                    const Zj = @as(T, @floatFromInt(self.atoms[j]));

                    const xj = self.coors[3 * j + 0];
                    const yj = self.coors[3 * j + 1];
                    const zj = self.coors[3 * j + 2];

                    const dx = xi - xj;
                    const dy = yi - yj;
                    const dz = zi - zj;

                    energy += (Zi * Zj) / std.math.sqrt(dx * dx + dy * dy + dz * dz);
                }
            }

            return energy;
        }

        /// Computes the one-electron nuclear-electron attraction potential matrix in the molecular basis.
        pub fn nuclear(self: @This(), gpa: Allocator) !Matrix(T) {
            const I = try Matrix(T).initZero(self.nbf, self.nbf, gpa);
            errdefer I.deinit(gpa);

            libint.libint_nuclear(I.data.ptr, self.ptr);

            return I;
        }

        /// Computes the first-order derivative of the nuclear attraction matrix with respect to nuclear coordinates.
        pub fn nuclearD1(self: @This(), gpa: Allocator) !Tensor(T, 3) {
            const shape = .{ 3 * libint.libint_nat(self.ptr), self.nbf, self.nbf };

            const I = try Tensor(T, 3).initZero(shape, gpa);
            errdefer I.deinit(gpa);

            libint.libint_nuclear_deriv(I.data.ptr, self.ptr);

            return I;
        }

        /// Computes the spin-blocked first-order derivative of the nuclear attraction matrix.
        pub fn nuclearD1Spin(self: @This(), gpa: Allocator) !Tensor(T, 3) {
            var V = try self.nuclearD1(gpa);
            defer V.deinit(gpa);

            const shape = .{ 3 * libint.libint_nat(self.ptr), 2 * self.nbf, 2 * self.nbf };

            var I = try Tensor(T, 3).initZero(shape, gpa);
            errdefer I.deinit(gpa);

            for (0..I.shape[0]) |k| for (0..self.nbf) |i| for (0..self.nbf) |j| {
                I.ptr(.{ k, i, j }).* = V.at(.{ k, i, j });

                I.ptr(.{ k, i + self.nbf, j + self.nbf }).* = I.at(.{ k, i, j });
            };

            return I;
        }

        /// Constructs the spin-blocked nuclear-electron attraction potential matrix.
        pub fn nuclearSpin(self: @This(), gpa: Allocator) !Matrix(T) {
            var V = try self.nuclear(gpa);
            defer V.deinit(gpa);

            var I = try Matrix(T).initZero(2 * self.nbf, 2 * self.nbf, gpa);
            errdefer I.deinit(gpa);

            for (0..self.nbf) |i| for (0..self.nbf) |j| {
                I.ptr(i, j).* = V.at(i, j);

                I.ptr(i + self.nbf, j + self.nbf).* = I.at(i, j);
            };

            return I;
        }

        /// Computes the overlap matrix elements representing the non-orthogonality of the spatial basis functions.
        pub fn overlap(self: @This(), gpa: Allocator) !Matrix(T) {
            const I = try Matrix(T).initZero(self.nbf, self.nbf, gpa);
            errdefer I.deinit(gpa);

            libint.libint_overlap(I.data.ptr, self.ptr);

            return I;
        }

        /// Computes the first-order derivative of the basis function overlap matrix with respect to nuclear coordinates.
        pub fn overlapD1(self: @This(), gpa: Allocator) !Tensor(T, 3) {
            const shape = .{ 3 * libint.libint_nat(self.ptr), self.nbf, self.nbf };

            const I = try Tensor(T, 3).initZero(shape, gpa);
            errdefer I.deinit(gpa);

            libint.libint_overlap_deriv(I.data.ptr, self.ptr);

            return I;
        }

        /// Computes the spin-blocked first-order derivative of the basis function overlap matrix.
        pub fn overlapD1Spin(self: @This(), gpa: Allocator) !Tensor(T, 3) {
            var S = try self.overlapD1(gpa);
            defer S.deinit(gpa);

            const shape = .{ 3 * libint.libint_nat(self.ptr), 2 * self.nbf, 2 * self.nbf };

            var I = try Tensor(T, 3).initZero(shape, gpa);
            errdefer I.deinit(gpa);

            for (0..I.shape[0]) |k| for (0..self.nbf) |i| for (0..self.nbf) |j| {
                I.ptr(.{ k, i, j }).* = S.at(.{ k, i, j });

                I.ptr(.{ k, i + self.nbf, j + self.nbf }).* = I.at(.{ k, i, j });
            };

            return I;
        }

        /// Constructs the spin-blocked basis function overlap matrix.
        pub fn overlapSpin(self: @This(), gpa: Allocator) !Matrix(T) {
            var S = try self.overlap(gpa);
            defer S.deinit(gpa);

            var I = try Matrix(T).initZero(2 * self.nbf, 2 * self.nbf, gpa);
            errdefer I.deinit(gpa);

            for (0..self.nbf) |i| for (0..self.nbf) |j| {
                I.ptr(i, j).* = S.at(i, j);

                I.ptr(i + self.nbf, j + self.nbf).* = I.at(i, j);
            };

            return I;
        }
    };
}
