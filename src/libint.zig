const std = @import("std");
const libint = @import("libint");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;

pub fn MolecularSystem(comptime T: type) type {
    return struct {
        ptr: *libint.SystemData,

        nbf: usize,

        pub fn init(system_path: [:0]const u8, basis_name: [:0]const u8) !@This() {
            const ptr = libint.init(system_path.ptr, basis_name.ptr) orelse {
                return error.InitializationFailed;
            };

            return .{ .ptr = ptr, .nbf = libint.nbf(ptr) };
        }

        pub fn deinit(self: *@This()) void {
            libint.deinit(self.ptr);
        }

        pub fn overlap(self: @This(), gpa: Allocator) !Matrix(T) {
            const I = try Matrix(T).initZero(self.nbf, self.nbf, gpa);

            libint.overlap(I.data.ptr, self.ptr);

            return I;
        }

        pub fn kinetic(self: @This(), gpa: Allocator) !Matrix(T) {
            const I = try Matrix(T).initZero(self.nbf, self.nbf, gpa);

            libint.kinetic(I.data.ptr, self.ptr);

            return I;
        }

        pub fn nuclear(self: @This(), gpa: Allocator) !Matrix(T) {
            const I = try Matrix(T).initZero(self.nbf, self.nbf, gpa);

            libint.nuclear(I.data.ptr, self.ptr);

            return I;
        }

        pub fn coulomb(self: @This(), gpa: Allocator) !Matrix(T) {
            const I = try Matrix(T).initZero(self.nbf * self.nbf, self.nbf * self.nbf, gpa);

            libint.coulomb(I.data.ptr, self.ptr);

            return I;
        }
    };
}
