const std = @import("std");

const libxc = @import("cimport.zig").libxc;

const primType = @import("value.zig").primType;

pub fn XCInput(comptime T: type) type {
    return struct {
        rho_val: []const T,

        sig_val: ?[]const T = null,
        tau_val: ?[]const T = null,
        lap_val: ?[]const T = null,
    };
}

pub fn XCOutput(comptime T: type) type {
    return struct {
        exc: []T,

        rho_pot: []T,

        sig_pot: ?[]T = null,
        tau_pot: ?[]T = null,
        lap_pot: ?[]T = null,
    };
}

pub fn evaluateXCFunctional(comptime T: type, out: XCOutput(T), func: libxc.xc_func_type, polarized: bool, inp: XCInput(T)) void {
    if (primType(T) != f64) @compileError("FUNCTIONAL EVALUATION NOW ONLY SUPPORTS F64 NUMBERS");

    if (func.info.*.family == libxc.XC_FAMILY_LDA or func.info.*.family == libxc.XC_FAMILY_HYB_LDA) {
        return evaluateXCLDA(T, out, func, polarized, inp);
    }

    if (func.info.*.family == libxc.XC_FAMILY_GGA or func.info.*.family == libxc.XC_FAMILY_HYB_GGA) {
        return evaluateXCGGA(T, out, func, polarized, inp);
    }

    if (func.info.*.family == libxc.XC_FAMILY_MGGA or func.info.*.family == libxc.XC_FAMILY_HYB_MGGA) {
        return evaluateXCMGGA(T, out, func, polarized, inp);
    }

    unreachable;
}

pub fn evaluateXCGGA(comptime T: type, out: XCOutput(T), func: libxc.xc_func_type, polarized: bool, inp: XCInput(T)) void {
    const CHUNK_SIZE = 512;

    const sig_val_unwrapped = inp.sig_val.?;
    const sig_pot_unwrapped = out.sig_pot.?;

    var e_xc_tmp: [1 * CHUNK_SIZE]T = undefined;
    var vrho_tmp: [2 * CHUNK_SIZE]T = undefined;
    var vsig_tmp: [3 * CHUNK_SIZE]T = undefined;

    var i: usize = 0;

    while (i < out.exc.len) : (i += @min(CHUNK_SIZE, out.exc.len - i)) {
        const batch_len = @min(CHUNK_SIZE, out.exc.len - i);

        const vrho_factor: usize = if (polarized) 2 else 1;
        const vsig_factor: usize = if (polarized) 3 else 1;

        const e_xc_size, const vrho_size = .{ batch_len, batch_len * vrho_factor };

        const e_xc_batch = e_xc_tmp[0..e_xc_size];
        const vrho_batch = vrho_tmp[0..vrho_size];

        @memset(e_xc_batch, 0);
        @memset(vrho_batch, 0);

        const vsig_batch = vsig_tmp[0 .. batch_len * vsig_factor];

        @memset(vsig_batch, 0);

        if (polarized) {
            const p1, const p2 = .{ inp.rho_val[2 * i ..].ptr, sig_val_unwrapped[3 * i ..].ptr };

            const p3 = e_xc_batch.ptr;
            const p4 = vrho_batch.ptr;
            const p5 = vsig_batch.ptr;

            libxc.xc_gga_exc_vxc(&func, batch_len, p1, p2, p3, p4, p5);

            for (0..batch_len) |j| {
                out.exc[i + j] += e_xc_batch[j];

                out.rho_pot[2 * (i + j) + 0] += vrho_batch[2 * j + 0];
                out.rho_pot[2 * (i + j) + 1] += vrho_batch[2 * j + 1];

                sig_pot_unwrapped[3 * (i + j) + 0] += vsig_batch[3 * j + 0];
                sig_pot_unwrapped[3 * (i + j) + 1] += vsig_batch[3 * j + 1];
                sig_pot_unwrapped[3 * (i + j) + 2] += vsig_batch[3 * j + 2];
            }
        }

        if (!polarized) {
            const p1, const p2 = .{ inp.rho_val[i..].ptr, sig_val_unwrapped[i..].ptr };

            const p3 = e_xc_batch.ptr;
            const p4 = vrho_batch.ptr;
            const p5 = vsig_batch.ptr;

            libxc.xc_gga_exc_vxc(&func, batch_len, p1, p2, p3, p4, p5);

            for (0..batch_len) |j| {
                out.exc[i + j], out.rho_pot[i + j] = .{ out.exc[i + j] + e_xc_batch[j], out.rho_pot[i + j] + vrho_batch[j] };

                sig_pot_unwrapped[i + j] += vsig_batch[j];
            }
        }
    }
}

pub fn evaluateXCLDA(comptime T: type, out: XCOutput(T), func: libxc.xc_func_type, polarized: bool, inp: XCInput(T)) void {
    const CHUNK_SIZE = 512;

    var e_xc_tmp: [1 * CHUNK_SIZE]T = undefined;
    var vrho_tmp: [2 * CHUNK_SIZE]T = undefined;

    var i: usize = 0;

    while (i < out.exc.len) : (i += @min(CHUNK_SIZE, out.exc.len - i)) {
        const batch_len = @min(CHUNK_SIZE, out.exc.len - i);

        const exc_size, const vrh_size = .{ batch_len, batch_len * (if (polarized) @as(usize, 2) else 1) };

        const exc_batch = e_xc_tmp[0..exc_size];
        const vrh_batch = vrho_tmp[0..vrh_size];

        @memset(exc_batch, 0);
        @memset(vrh_batch, 0);

        if (polarized) {
            libxc.xc_lda_exc_vxc(&func, batch_len, inp.rho_val[2 * i ..].ptr, exc_batch.ptr, vrh_batch.ptr);

            for (0..batch_len) |j| {
                out.exc[i + j] += exc_batch[j];

                out.rho_pot[2 * (i + j) + 0] += vrh_batch[2 * j + 0];
                out.rho_pot[2 * (i + j) + 1] += vrh_batch[2 * j + 1];
            }
        }

        if (!polarized) {
            libxc.xc_lda_exc_vxc(&func, batch_len, inp.rho_val[i..].ptr, exc_batch.ptr, vrh_batch.ptr);

            for (0..batch_len) |j| {
                out.exc[i + j] += exc_batch[j];

                out.rho_pot[i + j] += vrh_batch[j];
            }
        }
    }
}

pub fn evaluateXCMGGA(comptime T: type, out: XCOutput(T), func: libxc.xc_func_type, polarized: bool, inp: XCInput(T)) void {
    const CHUNK_SIZE = 512;

    const sig_val_unwrapped = inp.sig_val.?;
    const sig_pot_unwrapped = out.sig_pot.?;
    const tau_val_unwrapped = inp.tau_val.?;
    const tau_pot_unwrapped = out.tau_pot.?;
    const lap_val_unwrapped = inp.lap_val.?;
    const lap_pot_unwrapped = out.lap_pot.?;

    var e_xc_tmp: [1 * CHUNK_SIZE]T = undefined;
    var vrho_tmp: [2 * CHUNK_SIZE]T = undefined;
    var vsig_tmp: [3 * CHUNK_SIZE]T = undefined;
    var vtau_tmp: [2 * CHUNK_SIZE]T = undefined;
    var lapl_tmp: [2 * CHUNK_SIZE]T = undefined;
    var vlap_tmp: [2 * CHUNK_SIZE]T = undefined;

    var i: usize = 0;

    while (i < out.exc.len) : (i += @min(CHUNK_SIZE, out.exc.len - i)) {
        const batch_len = @min(CHUNK_SIZE, out.exc.len - i);

        const vrho_factor: usize = if (polarized) 2 else 1;
        const vtau_factor: usize = if (polarized) 2 else 1;
        const vsig_factor: usize = if (polarized) 3 else 1;

        const e_xc_size, const vrho_size = .{ batch_len, batch_len * vrho_factor };

        const vsig_size = batch_len * vsig_factor;
        const vtau_size = batch_len * vtau_factor;

        const e_xc_batch = e_xc_tmp[0..e_xc_size];
        const vrho_batch = vrho_tmp[0..vrho_size];

        @memset(e_xc_batch, 0);
        @memset(vrho_batch, 0);

        const vsig_batch = vsig_tmp[0..vsig_size];
        const vtau_batch = vtau_tmp[0..vtau_size];

        @memset(vsig_batch, 0);
        @memset(vtau_batch, 0);

        const lapl_batch = lapl_tmp[0..vrho_size];
        const vlap_batch = vlap_tmp[0..vrho_size];

        @memset(vlap_batch, 0);

        if (polarized) {
            const p1 = inp.rho_val[2 * i ..].ptr;

            const p2 = sig_val_unwrapped[3 * i ..].ptr;
            const p4 = tau_val_unwrapped[2 * i ..].ptr;

            @memcpy(lapl_batch, lap_val_unwrapped[2 * i .. 2 * i + vrho_size]);

            const p3 = lapl_batch.ptr;
            const p5 = e_xc_batch.ptr;
            const p6 = vrho_batch.ptr;
            const p7 = vsig_batch.ptr;
            const p8 = vlap_batch.ptr;
            const p9 = vtau_batch.ptr;

            libxc.xc_mgga_exc_vxc(&func, batch_len, p1, p2, p3, p4, p5, p6, p7, p8, p9);

            for (0..batch_len) |j| {
                out.exc[i + j] += e_xc_batch[j];

                out.rho_pot[2 * (i + j) + 0] += vrho_batch[2 * j + 0];
                out.rho_pot[2 * (i + j) + 1] += vrho_batch[2 * j + 1];

                sig_pot_unwrapped[3 * (i + j) + 0] += vsig_batch[3 * j + 0];
                sig_pot_unwrapped[3 * (i + j) + 1] += vsig_batch[3 * j + 1];
                sig_pot_unwrapped[3 * (i + j) + 2] += vsig_batch[3 * j + 2];
                tau_pot_unwrapped[2 * (i + j) + 0] += vtau_batch[2 * j + 0];
                tau_pot_unwrapped[2 * (i + j) + 1] += vtau_batch[2 * j + 1];
                lap_pot_unwrapped[2 * (i + j) + 0] += vlap_batch[2 * j + 0];
                lap_pot_unwrapped[2 * (i + j) + 1] += vlap_batch[2 * j + 1];
            }
        }

        if (!polarized) {
            const p1 = inp.rho_val[i..].ptr;

            const p2 = sig_val_unwrapped[i..].ptr;
            const p4 = tau_val_unwrapped[i..].ptr;

            @memcpy(lapl_batch, lap_val_unwrapped[i .. i + vrho_size]);

            const p3 = lapl_batch.ptr;
            const p5 = e_xc_batch.ptr;
            const p6 = vrho_batch.ptr;
            const p7 = vsig_batch.ptr;
            const p8 = vlap_batch.ptr;
            const p9 = vtau_batch.ptr;

            libxc.xc_mgga_exc_vxc(&func, batch_len, p1, p2, p3, p4, p5, p6, p7, p8, p9);

            for (0..batch_len) |j| {
                out.exc[i + j], out.rho_pot[i + j] = .{ out.exc[i + j] + e_xc_batch[j], out.rho_pot[i + j] + vrho_batch[j] };

                sig_pot_unwrapped[i + j] += vsig_batch[j];
                tau_pot_unwrapped[i + j] += vtau_batch[j];
                lap_pot_unwrapped[i + j] += vlap_batch[j];
            }
        }
    }
}
