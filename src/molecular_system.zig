const std = @import("std");

const libint = @import("cimport.zig").libint;

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

pub fn MolecularSystem(comptime T: type) type {
    return struct {
        ptr: *libint.SystemData,

        atoms: []i32,
        coors: []f64,
        bf2at: []i32,

        nbf: usize,
        nel: usize,

        pub fn init(system: []const u8, basis: []const u8, gpa: Allocator) !@This() {
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

            const coors = try gpa.alloc(f64, 3 * nat);
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

            return .{ .ptr = ptr, .nbf = nbf, .atoms = atoms, .coors = coors, .bf2at = bf2at, .nel = nel };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            gpa.free(self.atoms);
            gpa.free(self.coors);
            gpa.free(self.bf2at);

            libint.libint_deinit(self.ptr);
        }

        pub fn coulomb(self: @This(), gpa: Allocator) !Tensor(T, 4) {
            const I = try Tensor(T, 4).initZero(.{ self.nbf, self.nbf, self.nbf, self.nbf }, gpa);
            errdefer I.deinit(gpa);

            libint.libint_coulomb(I.data.ptr, self.ptr);

            return I;
        }

        pub fn coulombD1(self: @This(), gpa: Allocator) !Tensor(T, 5) {
            const shape = .{ 3 * libint.libint_nat(self.ptr), self.nbf, self.nbf, self.nbf, self.nbf };

            const I = try Tensor(T, 5).initZero(shape, gpa);
            errdefer I.deinit(gpa);

            libint.libint_coulomb_deriv(I.data.ptr, self.ptr);

            return I;
        }

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

        pub fn kinetic(self: @This(), gpa: Allocator) !Matrix(T) {
            const I = try Matrix(T).initZero(self.nbf, self.nbf, gpa);
            errdefer I.deinit(gpa);

            libint.libint_kinetic(I.data.ptr, self.ptr);

            return I;
        }

        pub fn kineticD1(self: @This(), gpa: Allocator) !Tensor(T, 3) {
            const shape = .{ 3 * libint.libint_nat(self.ptr), self.nbf, self.nbf };

            const I = try Tensor(T, 3).initZero(shape, gpa);
            errdefer I.deinit(gpa);

            libint.libint_kinetic_deriv(I.data.ptr, self.ptr);

            return I;
        }

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

        pub fn nuclear(self: @This(), gpa: Allocator) !Matrix(T) {
            const I = try Matrix(T).initZero(self.nbf, self.nbf, gpa);
            errdefer I.deinit(gpa);

            libint.libint_nuclear(I.data.ptr, self.ptr);

            return I;
        }

        pub fn nuclearD1(self: @This(), gpa: Allocator) !Tensor(T, 3) {
            const shape = .{ 3 * libint.libint_nat(self.ptr), self.nbf, self.nbf };

            const I = try Tensor(T, 3).initZero(shape, gpa);
            errdefer I.deinit(gpa);

            libint.libint_nuclear_deriv(I.data.ptr, self.ptr);

            return I;
        }

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

        pub fn overlap(self: @This(), gpa: Allocator) !Matrix(T) {
            const I = try Matrix(T).initZero(self.nbf, self.nbf, gpa);
            errdefer I.deinit(gpa);

            libint.libint_overlap(I.data.ptr, self.ptr);

            return I;
        }

        pub fn overlapD1(self: @This(), gpa: Allocator) !Tensor(T, 3) {
            const shape = .{ 3 * libint.libint_nat(self.ptr), self.nbf, self.nbf };

            const I = try Tensor(T, 3).initZero(shape, gpa);
            errdefer I.deinit(gpa);

            libint.libint_overlap_deriv(I.data.ptr, self.ptr);

            return I;
        }

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
