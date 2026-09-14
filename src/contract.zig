//! Tensor contraction engine utilizing BLAS matrix multiplication for multi-dimensional physics arrays.

const std = @import("std");

const mmSlice = @import("linear_algebra.zig").mmSlice;
const primType = @import("value.zig").primType;

/// Specification plan storing tensor index permutations and contraction dimensionality.
fn ContractionPlan(comptime spec: []const u8) type {
    const comma = std.mem.indexOfScalar(u8, spec, ',') orelse @compileError("MISSING ',' IN SPEC");
    const arrow = std.mem.indexOf(u8, spec, "->") orelse @compileError("MISSING '->' SYM IN SPEC");

    const dim_a = comma;

    const dim_b, const dim_c = .{ arrow - comma - 1, spec.len - arrow - 2 };

    return struct {
        pub const na = dim_a;
        pub const nb = dim_b;
        pub const nc = dim_c;

        perm_a: [na]usize,
        perm_b: [nb]usize,
        perm_c: [nc]usize,

        num_k: usize,
    };
}

/// Contracts two tensors according to an Einstein summation specification using CBLAS dgemm.
pub fn contract(comptime T: type, comptime spec: []const u8, C: anytype, A: anytype, B: anytype, alpha: T, beta: T, gpa: ?std.mem.Allocator) !void {
    if (comptime primType(T) != f64) @compileError("CONTRACT ONLY SUPPORTS F64 NUMBERS");

    const a_ten, const b_ten, const c_ten = .{ A.asTensor(), B.asTensor(), C.asTensor() };

    const plan = comptime parseContraction(spec);

    const na, const nb, const nc = .{ @TypeOf(plan).na, @TypeOf(plan).nb, @TypeOf(plan).nc };

    std.debug.assert(a_ten.shape.len == na and b_ten.shape.len == nb and c_ten.shape.len == nc);

    const num_i, const num_j = .{ na - plan.num_k, nb - plan.num_k };

    var a_shape: [na]usize, var b_shape: [nb]usize = .{ undefined, undefined };

    var n: usize, var m: usize, var k: usize = .{ 1, 1, 1 };

    inline for (0..na) |i| {
        a_shape[i] = a_ten.shape[plan.perm_a[i]];

        if (i < num_i) m *= a_shape[i] else k *= a_shape[i];
    }

    inline for (0..nb) |i| {
        b_shape[i] = b_ten.shape[plan.perm_b[i]];

        if (i >= plan.num_k) n *= b_shape[i] else std.debug.assert(a_shape[num_i + i] == b_shape[i]);
    }

    var shape_gemm: [nc]usize = undefined;

    inline for (0..num_i) |i| shape_gemm[i] = a_shape[i];

    inline for (0..num_j) |i| {
        shape_gemm[num_i + i] = b_shape[plan.num_k + i];
    }

    inline for (0..nc) |i| std.debug.assert(c_ten.shape[i] == shape_gemm[plan.perm_c[i]]);

    const a_slice, const a_buf = try prepareSlice(T, na, a_ten, plan.perm_a, gpa);
    defer if (a_buf) |buf| (gpa.?).free(buf);

    const b_slice, const b_buf = try prepareSlice(T, nb, b_ten, plan.perm_b, gpa);
    defer if (b_buf) |buf| (gpa.?).free(buf);

    if (isIdentity(plan.perm_c)) {
        return mmSlice(T, c_ten.data, a_slice, b_slice, m, n, k, alpha, beta, false, false);
    }

    const alloc = gpa orelse return error.AllocatorRequired;

    const c_buf = try alloc.alloc(T, c_ten.data.len);
    defer alloc.free(c_buf);

    mmSlice(T, c_buf, a_slice, b_slice, m, n, k, alpha, 0.0, false, false);

    if (beta == 0) {
        return permuteSlice(T, nc, c_ten.data, c_buf, shape_gemm, plan.perm_c);
    }

    const c_perm = try alloc.alloc(T, c_ten.data.len);
    defer alloc.free(c_perm);

    permuteSlice(T, nc, c_perm, c_buf, shape_gemm, plan.perm_c);

    for (c_ten.data, c_perm) |*c, p| c.* = p + beta * c.*;
}

/// Checks if an index permutation vector is an identity mapping.
fn isIdentity(comptime order: anytype) bool {
    inline for (order, 0..) |val, i| {
        if (val != i) return false;
    }

    return true;
}

/// Parses an Einstein summation specification into index permutations and contraction count.
fn parseContraction(comptime spec: []const u8) ContractionPlan(spec) {
    const Plan = ContractionPlan(spec);

    const na, const nb, var num_k = .{ Plan.na, Plan.nb, 0 };

    const idx_a, const idx_b, const idx_c = .{ spec[0..na], spec[na + 1 .. na + 1 + nb], spec[na + nb + 3 ..] };

    for (idx_a) |ca| {
        if (std.mem.indexOfScalar(u8, idx_b, ca) != null) num_k += 1;
    }

    const num_i = Plan.na - num_k;
    const num_j = Plan.nb - num_k;

    if (num_i + num_j != Plan.nc) @compileError("FREE INDICES COUNT MUST MATCH RANK OF OUTPUT TENSOR");

    var perm_a: [Plan.na]usize, var free_a: usize, var cont_a: usize = .{ undefined, 0, num_i };

    for (idx_a, 0..) |ca, i| {
        const has_cont = std.mem.indexOfScalar(u8, idx_b, ca) != null;

        if (has_cont) {
            perm_a[cont_a], cont_a = .{ i, cont_a + 1 };
        }

        if (!has_cont) {
            perm_a[free_a], free_a = .{ i, free_a + 1 };
        }
    }

    var perm_b: [Plan.nb]usize = undefined;

    for (0..num_k) |k| {
        perm_b[k] = std.mem.indexOfScalar(u8, idx_b, idx_a[perm_a[num_i + k]]).?;
    }

    var free_b: usize = num_k;

    for (idx_b, 0..) |cb, j| if (std.mem.indexOfScalar(u8, idx_a, cb) == null) {
        perm_b[free_b], free_b = .{ j, free_b + 1 };
    };

    var gemm_idx: [Plan.nc]u8 = undefined;

    for (0..num_i) |i| {
        gemm_idx[i] = idx_a[perm_a[i]];
    }

    for (0..num_j) |j| {
        gemm_idx[num_i + j] = idx_b[perm_b[num_k + j]];
    }

    var perm_c: [Plan.nc]usize = undefined;

    for (idx_c, 0..) |cc, i| {
        perm_c[i] = std.mem.indexOfScalar(u8, &gemm_idx, cc) orelse @compileError("INDEX IN OUTPUT NOT FOUND IN FREE INDICES");
    }

    return .{
        .perm_a = perm_a,
        .perm_b = perm_b,
        .perm_c = perm_c,

        .num_k = num_k,
    };
}

/// Permutes multi-dimensional tensor data in-place or out-of-place according to specified axis ordering.
fn permuteSlice(comptime T: type, comptime N: usize, dst: []T, src: []const T, shape: [N]usize, comptime order: [N]usize) void {
    if (comptime N == 0 or isIdentity(order)) {
        @memcpy(dst, src);

        return;
    }

    var strides: [N]usize, var str: usize = .{ undefined, 1 };

    inline for (0..N) |k| {
        const j = N - 1 - k;

        strides[order[j]], str = .{ str, str * shape[order[j]] };
    }

    var idx, var offset: usize = .{ std.mem.zeroes([N]usize), 0 };

    for (src) |val| {
        dst[offset] = val;

        inline for (0..N) |k| {
            const dim = N - 1 - k;

            idx[dim] += 1;

            if (idx[dim] < shape[dim]) {
                offset += strides[dim];

                break;
            }

            if (idx[dim] == shape[dim]) {
                idx[dim], offset = .{ 0, offset - strides[dim] * (shape[dim] - 1) };
            }
        }
    }
}

/// Prepares a contiguous slice for GEMM by permuting tensor data if axes are not in identity order.
fn prepareSlice(comptime T: type, comptime N: usize, ten: anytype, comptime perm: [N]usize, gpa: ?std.mem.Allocator) !struct { []const T, ?[]T } {
    if (isIdentity(perm)) return .{ ten.data, null };

    const alloc = gpa orelse return error.AllocatorRequired;

    const buf = try alloc.alloc(T, ten.data.len);

    permuteSlice(T, N, buf, ten.data, ten.shape, perm);

    return .{ buf, buf };
}
