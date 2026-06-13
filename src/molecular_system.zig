const libint = @import("libint");
const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

pub fn MolecularSystem(comptime T: type) type {
    return struct {
        ptr: *libint.SystemData,

        atoms: []i32,
        coors: []f64,

        nbf: usize,
        nel: usize,

        pub fn init(system: []const u8, basis: []const u8, gpa: Allocator) !@This() {
            const sys_c = try gpa.dupeSentinel(u8, system, 0);
            defer gpa.free(sys_c);

            const bas_c = try gpa.dupeSentinel(u8, basis, 0);
            defer gpa.free(bas_c);

            const ptr = libint.init(sys_c.ptr, bas_c.ptr) orelse {
                return error.InitializationFailed;
            };

            const nat = libint.nat(ptr);

            const atoms = try gpa.alloc(i32, 1 * nat);
            const coors = try gpa.alloc(f64, 3 * nat);

            libint.atoms(atoms.ptr, ptr);
            libint.coors(coors.ptr, ptr);

            var nel: usize = 0;

            for (atoms) |z| {
                nel += @intCast(z);
            }

            return .{ .ptr = ptr, .nbf = libint.nbf(ptr), .atoms = atoms, .coors = coors, .nel = nel };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            gpa.free(self.atoms);
            gpa.free(self.coors);

            libint.deinit(self.ptr);
        }

        pub fn overlap(self: @This(), gpa: Allocator) !Matrix(T) {
            const I = try Matrix(T).initZero(self.nbf, self.nbf, gpa);
            errdefer I.deinit(gpa);

            libint.overlap(I.data.ptr, self.ptr);

            return I;
        }

        pub fn kinetic(self: @This(), gpa: Allocator) !Matrix(T) {
            const I = try Matrix(T).initZero(self.nbf, self.nbf, gpa);
            errdefer I.deinit(gpa);

            libint.kinetic(I.data.ptr, self.ptr);

            return I;
        }

        pub fn nuclear(self: @This(), gpa: Allocator) !Matrix(T) {
            const I = try Matrix(T).initZero(self.nbf, self.nbf, gpa);
            errdefer I.deinit(gpa);

            libint.nuclear(I.data.ptr, self.ptr);

            return I;
        }

        pub fn coulomb(self: @This(), gpa: Allocator) !Tensor(T, 4) {
            const I = try Tensor(T, 4).initZero(.{ self.nbf, self.nbf, self.nbf, self.nbf }, gpa);
            errdefer I.deinit(gpa);

            libint.coulomb(I.data.ptr, self.ptr);

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

        pub fn coulombSpin(self: @This(), gpa: Allocator) !Tensor(T, 4) {
            var J = try self.coulomb(gpa);
            defer J.deinit(gpa);

            var I = try Tensor(T, 4).initZero(.{ 2 * self.nbf, 2 * self.nbf, 2 * self.nbf, 2 * self.nbf }, gpa);
            errdefer I.deinit(gpa);

            for (0..self.nbf) |ma| for (0..self.nbf) |mb| for (0..self.nbf) |mc| for (0..self.nbf) |md| {
                const val = J.at(.{ ma, mb, mc, md });

                I.ptr(.{ ma, mb, mc, md }).* = val;

                I.ptr(.{ ma + self.nbf, mb, mc + self.nbf, md }).* = val;
                I.ptr(.{ ma, mb + self.nbf, mc, md + self.nbf }).* = val;

                I.ptr(.{ ma + self.nbf, mb + self.nbf, mc + self.nbf, md + self.nbf }).* = val;
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
    };
}
