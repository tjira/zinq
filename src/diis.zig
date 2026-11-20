//! Implements Direct Inversion in the Iterative Subspace (DIIS) acceleration for SCF and response equations.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Vector = @import("tensor.zig").Vector;

const axpy = @import("linear_algebra.zig").axpy;
const dot = @import("linear_algebra.zig").dot;
const getrf = @import("linear_algebra.zig").getrf;
const getrs = @import("linear_algebra.zig").getrs;

/// Manages iterative subspace histories and extrapolates solution vectors via DIIS error minimization.
pub fn Diis(comptime T: type) type {
    return struct {
        fck_hist: std.ArrayList(Vector(T)),
        err_hist: std.ArrayList(Vector(T)),

        max_diis: usize,

        /// Initializes the DIIS accelerator with a fixed history subspace capacity.
        pub fn init(max_diis: usize) @This() {
            return .{ .fck_hist = .empty, .err_hist = .empty, .max_diis = max_diis };
        }

        /// Frees all history vectors and deallocates the subspace storage.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            for (self.fck_hist.items) |*v| v.deinit(gpa);
            for (self.err_hist.items) |*v| v.deinit(gpa);

            self.fck_hist.deinit(gpa);
            self.err_hist.deinit(gpa);
        }

        /// Clears the subspace histories and resets the accelerator state upon numerical failure.
        pub fn reset(self: *@This(), gpa: Allocator) void {
            for (self.fck_hist.items) |*v| v.deinit(gpa);
            for (self.err_hist.items) |*v| v.deinit(gpa);

            self.fck_hist.clearRetainingCapacity();
            self.err_hist.clearRetainingCapacity();
        }

        /// Appends an iterate and error vector, updating the target via subspace extrapolation.
        pub fn update(self: *@This(), target: anytype, err: anytype, symmetric: bool, gpa: Allocator) !void {
            var target_vec, var err_vec = .{ Vector(T).fromSlice(target.data), Vector(T).fromSlice(err.data) };

            std.debug.assert(target_vec.length() == err_vec.length());

            try self.fck_hist.ensureUnusedCapacity(gpa, 1);
            try self.err_hist.ensureUnusedCapacity(gpa, 1);

            var f_copy = try Vector(T).init(target_vec.length(), gpa);
            errdefer f_copy.deinit(gpa);

            @memcpy(f_copy.data, target_vec.data);

            if (self.fck_hist.items.len >= self.max_diis) {
                var old_f = self.fck_hist.orderedRemove(0);
                var old_e = self.err_hist.orderedRemove(0);

                old_f.deinit(gpa);
                old_e.deinit(gpa);
            }

            self.fck_hist.appendAssumeCapacity(f_copy);
            self.err_hist.appendAssumeCapacity(err_vec);

            if (self.fck_hist.items.len >= 2) {
                extrapolate(T, self.fck_hist.items, self.err_hist.items, &target_vec, gpa) catch {
                    const last_f = self.fck_hist.pop().?;
                    const last_e = self.err_hist.pop().?;

                    self.reset(gpa);

                    self.fck_hist.appendAssumeCapacity(last_f);
                    self.err_hist.appendAssumeCapacity(last_e);
                };
            }

            if (symmetric) {
                if (comptime @hasDecl(@TypeOf(target.*), "symmetrize")) {
                    target.symmetrize();
                }
            }
        }

        /// Computes linear combination coefficients from the B matrix and extrapolates the target.
        fn extrapolate(comptime TF: type, fck_hist: []const Vector(TF), err_hist: []const Vector(TF), target_vec: *Vector(TF), gpa: Allocator) !void {
            const n = fck_hist.len;

            var B = try Matrix(TF).initZero(n + 1, n + 1, gpa);
            defer B.deinit(gpa);

            var b = try Matrix(TF).initZero(n + 1, 1, gpa);
            defer b.deinit(gpa);

            for (0..n) |i| for (i..n) |j| {
                const s = dot(TF, err_hist[i], err_hist[j]);

                B.ptr(i, j).* = s;
                B.ptr(j, i).* = s;
            };

            for (0..n) |i| {
                B.ptr(i, n).* = -1.0;
                B.ptr(n, i).* = -1.0;
            }

            b.ptr(n, 0).* = -1.0;

            const ipiv = try gpa.alloc(i32, n + 1);
            defer gpa.free(ipiv);

            try getrf(TF, &B, ipiv);

            var c = try Matrix(TF).init(n + 1, 1, gpa);
            defer c.deinit(gpa);

            try getrs(TF, &c, B, ipiv, b);

            target_vec.zero();

            for (0..n) |i| {
                axpy(TF, target_vec, fck_hist[i], c.at(i, 0));
            }
        }
    };
}
