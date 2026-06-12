const libint = @import("libint");
const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

pub fn MolecularSystem(comptime T: type) type {
    return struct {
        ptr: *libint.SystemData,

        nbf: usize,

        pub fn init(system: []const u8, basis: []const u8, gpa: Allocator) !@This() {
            const sys_c = try gpa.dupeSentinel(u8, system, 0);
            defer gpa.free(sys_c);

            const bas_c = try gpa.dupeSentinel(u8, basis, 0);
            defer gpa.free(bas_c);

            const ptr = libint.init(sys_c.ptr, bas_c.ptr) orelse {
                return error.InitializationFailed;
            };

            return .{ .ptr = ptr, .nbf = libint.nbf(ptr) };
        }

        pub fn deinit(self: *@This()) void {
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
    };
}
