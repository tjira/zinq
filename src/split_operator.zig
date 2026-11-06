//! Propagates multi-state quantum wavepackets using split-operator Fourier methods with absorbing boundaries.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Grid = @import("wavepacket.zig").Grid;
const Hamiltonian = @import("wavepacket.zig").Hamiltonian;
const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const Vector = @import("tensor.zig").Vector;

/// Wavepacket propagator using the split-operator Fourier method with absorbing boundary potentials.
pub fn SplitOperator(comptime T: type) type {
    return struct {
        R: ?Matrix(Complex(T)),
        K: ?Vector(Complex(T)),
        cap_weight: ?Vector(T),

        r_buf: ?[]Complex(T),

        dt: Complex(T),

        /// Initializes the kinetic and potential propagation operators and absorbing boundary weights.
        pub fn init(grid: Grid(T), ham: Hamiltonian(T), pot: Potential(T), capopt: anytype, dt: Complex(T), optimize_memory: bool, gpa: Allocator) !@This() {
            if (optimize_memory) {
                const r_buf = try gpa.alloc(Complex(T), pot.nstate() * pot.nstate());
                errdefer gpa.free(r_buf);

                return @This(){ .R = null, .K = null, .cap_weight = null, .r_buf = r_buf, .dt = dt };
            }

            const nstate = pot.nstate();

            var R = try Matrix(Complex(T)).init(grid.nrow(), nstate * nstate, gpa);
            errdefer R.deinit(gpa);

            var K = try Vector(Complex(T)).initZero(grid.nrow(), gpa);
            errdefer K.deinit(gpa);

            var cap_weight = try Vector(T).initZero(grid.nrow(), gpa);
            errdefer cap_weight.deinit(gpa);

            for (0..grid.nrow()) |i| {
                K.data[i] = std.math.complex.exp(Complex(T).init(0, -ham.getK(grid, i)).mul(dt));
            }

            var prop = @This(){ .R = R, .K = K, .cap_weight = cap_weight, .r_buf = null, .dt = dt };

            try prop.update(grid, ham, pot, 0, capopt);

            return prop;
        }

        /// Deallocates split-operator propagation matrices.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.r_buf) |buf| gpa.free(buf);

            if (self.R) |*R| R.deinit(gpa);
            if (self.K) |*K| K.deinit(gpa);

            if (self.cap_weight) |*cw| cw.deinit(gpa);
        }

        /// Retrieves the absorbing boundary decay coefficient value at a grid coordinate to damp outgoing flux.
        pub fn getCap(self: @This(), grid: Grid(T), capopt: anytype, i: usize) T {
            if (self.cap_weight) |cw| return cw.data[i];

            var cap_sum: T = 0;

            if (capopt) |cap| for (0..grid.ncol()) |j| {
                const r_val = grid.getR(i, j);

                const min_bound = cap.bounds[j][0];
                const max_bound = cap.bounds[j][1];

                const dist = @max(0, @max(min_bound - r_val, r_val - max_bound));

                if (dist > 0) {
                    cap_sum += std.math.exp(cap.exponent * dist) - 1;
                }
            };

            const cap_decay = std.math.exp(-0.5 * cap_sum * self.dt.magnitude());

            return 1 - cap_decay * cap_decay;
        }

        /// Retrieves the kinetic propagator phase factor representing momentum-dependent grid transitions.
        pub fn getK(self: @This(), grid: Grid(T), ham: Hamiltonian(T), i: usize) Complex(T) {
            if (self.K) |K| return K.at(i);

            return std.math.complex.exp(Complex(T).init(0, -ham.getK(grid, i)).mul(self.dt));
        }

        /// Retrieves the diabatic propagator row representing potential-induced transitions at a grid point.
        pub fn getR(self: @This(), grid: Grid(T), ham: Hamiltonian(T), pot: Potential(T), t: T, capopt: anytype, i: usize) ![]Complex(T) {
            if (self.R) |R| return R.rowSlice(i);

            const nstate, var cap_decay: T = .{ pot.nstate(), 1 };

            if (capopt != null) {
                cap_decay = std.math.sqrt(@max(0, 1 - self.getCap(grid, capopt, i)));
            }

            const w, const u, _ = try ham.getTriple(grid, pot, t, i);

            const buffer = self.r_buf.?;

            for (0..nstate) |k| for (0..nstate) |j| {
                var sum = Complex(T).init(0, 0);

                for (0..nstate) |m| {
                    const U_km = u[k * nstate + m];
                    const U_jm = u[j * nstate + m];

                    const phase = std.math.complex.exp(Complex(T).init(0, -0.5 * w[m]).mul(self.dt));

                    sum = sum.add(Complex(T).init(phase.re * U_km * U_jm, phase.im * U_km * U_jm));
                }

                buffer[k * nstate + j] = sum.mul(Complex(T).init(cap_decay, 0));
            };

            return buffer;
        }

        /// Propagates the wavepacket by one time step using potential and kinetic operator splits.
        pub fn step(self: @This(), sim: anytype, capopt: anytype, adia: bool, track_pop: bool, t: T, gpa: Allocator) !void {
            if (track_pop) try self.accumAbsorbed(sim, capopt, adia, t, gpa);

            try self.applyR(sim, capopt, t, gpa);

            try self.applyK(sim);

            if (track_pop) try self.accumAbsorbed(sim, capopt, adia, t, gpa);

            try self.applyR(sim, capopt, t, gpa);
        }

        /// Updates the potential propagator matrix and absorbing boundary exponential decay.
        pub fn update(self: *@This(), grid: Grid(T), ham: Hamiltonian(T), pot: Potential(T), t: T, capopt: anytype) !void {
            if (self.R == null) return;

            const nstate = pot.nstate();

            for (0..self.R.?.nrow()) |i| {
                var cap_sum: T = 0;

                if (capopt) |cap| for (0..grid.ncol()) |j| {
                    const r_val = grid.getR(i, j);

                    const min_bound = cap.bounds[j][0];
                    const max_bound = cap.bounds[j][1];

                    const dist = @max(0, @max(min_bound - r_val, r_val - max_bound));

                    if (dist > 0) {
                        cap_sum += std.math.exp(cap.exponent * dist) - 1;
                    }
                };

                const cap_decay = std.math.exp(-0.5 * cap_sum * self.dt.magnitude());

                self.cap_weight.?.ptr(i).* = 1 - cap_decay * cap_decay;

                const R_i = try self.getR(grid, ham, pot, t, capopt, i);

                const w, const u, _ = try ham.getTriple(grid, pot, t, i);

                for (0..nstate) |k| for (0..nstate) |j| {
                    var sum = Complex(T).init(0, 0);

                    for (0..nstate) |m| {
                        const U_km = u[k * nstate + m];
                        const U_jm = u[j * nstate + m];

                        const phase = std.math.complex.exp(Complex(T).init(0, -0.5 * w[m]).mul(self.dt));

                        sum = sum.add(Complex(T).init(phase.re * U_km * U_jm, phase.im * U_km * U_jm));
                    }

                    R_i[k * nstate + j] = sum.mul(Complex(T).init(cap_decay, 0));
                };
            }
        }

        /// Computes and accumulates the population absorbed by the boundary potential.
        fn accumAbsorbed(self: @This(), sim: anytype, capopt: anytype, adia: bool, t: T, gpa: Allocator) !void {
            const nstate = sim.wfn.W.nrow();

            const sums = try gpa.alloc(T, nstate);
            defer gpa.free(sums);

            @memset(sums, 0);

            for (0..sim.wfn.W.ncol()) |j| {
                const w = self.getCap(sim.wfn_kpgrids, capopt, j);

                if (w == 0) continue;

                if (adia) {
                    try accumAbsorbedAdia(T, sim, sums, w, t, j);
                }

                if (!adia) {
                    accumAbsorbedDia(T, sim, sums, w, j);
                }
            }

            for (0..nstate) |i| {
                sim.pop_apabs.ptr(i).* += sums[i] * sim.wfn_kpgrids.dr;
            }
        }

        /// Computes adiabatic absorbed boundary population contributions for a grid point.
        fn accumAbsorbedAdia(comptime U: type, sim: anytype, sums: []U, w: U, t: U, j: usize) !void {
            const nstate = sim.wfn.W.nrow();

            const triple = try sim.hams.getTriple(sim.wfn_kpgrids, sim.epoten, t, j);

            for (0..nstate) |i| {
                var adia_re: U = 0;
                var adia_im: U = 0;

                for (0..nstate) |k| {
                    const u = triple[1][k * nstate + i];

                    adia_re += sim.wfn.W.at(k, j).re * u;
                    adia_im += sim.wfn.W.at(k, j).im * u;
                }

                sums[i] += w * (adia_re * adia_re + adia_im * adia_im);
            }
        }

        /// Computes diabatic absorbed boundary population contributions for a grid point.
        fn accumAbsorbedDia(comptime U: type, sim: anytype, sums: []U, w: U, j: usize) void {
            const nstate = sim.wfn.W.nrow();

            for (0..nstate) |i| {
                sums[i] += w * sim.wfn.W.at(i, j).squaredMagnitude();
            }
        }

        /// Applies the kinetic energy propagator in momentum space using FFT.
        fn applyK(self: @This(), sim: anytype) !void {
            sim.wfn.fft(-1);

            defer {
                sim.wfn.fft(1);
            }

            for (0..sim.wfn.W.nrow()) |i| for (0..sim.wfn.W.ncol()) |j| {
                sim.wfn.W.ptr(i, j).* = self.getK(sim.wfn_kpgrids, sim.hams, j).mul(sim.wfn.W.at(i, j));
            };
        }

        /// Applies the potential energy propagator in position space.
        fn applyR(self: @This(), sim: anytype, capopt: anytype, t: T, gpa: Allocator) !void {
            const nstate = sim.wfn.W.nrow();

            if (nstate < 16) {
                inline for (1..16) |N| {
                    if (nstate == N) return try self.applyRFast(N, sim, capopt, t);
                }
            }

            try self.applyRSlow(sim, capopt, t, gpa);
        }

        /// Fast potential propagation in coordinate space for small state spaces using stack-allocated structures.
        fn applyRFast(self: @This(), comptime N: usize, sim: anytype, capopt: anytype, t: T) !void {
            return if (self.R) |_| self.applyRFastStatic(N, sim) else self.applyRFastDynamic(N, sim, capopt, t);
        }

        /// Fast potential propagation in coordinate space with on-the-fly propagator computation.
        fn applyRFastDynamic(self: @This(), comptime N: usize, sim: anytype, capopt: anytype, t: T) !void {
            var rows: [N][]Complex(T) = undefined;

            var w_col: [N]Complex(T) = undefined;

            inline for (0..N) |k| {
                rows[k] = sim.wfn.W.rowSlice(k);
            }

            const ncol = sim.wfn.W.ncol();

            for (0..ncol) |j| {
                const R_j = try self.getR(sim.wfn_kpgrids, sim.hams, sim.epoten, t, capopt, j);

                inline for (0..N) |k| {
                    w_col[k] = rows[k][j];
                }

                inline for (0..N) |i| {
                    var sum_re: T = 0;
                    var sum_im: T = 0;

                    inline for (0..N) |k| {
                        const r, const w = .{ R_j[i * N + k], w_col[k] };

                        sum_re += r.re * w.re - r.im * w.im;
                        sum_im += r.re * w.im + r.im * w.re;
                    }

                    rows[i][j] = Complex(T).init(sum_re, sum_im);
                }
            }
        }

        /// Fast potential propagation in coordinate space with precomputed static propagator matrices.
        fn applyRFastStatic(self: @This(), comptime N: usize, sim: anytype) void {
            var rows: [N][]Complex(T) = undefined;

            var w_col: [N]Complex(T) = undefined;

            inline for (0..N) |k| {
                rows[k] = sim.wfn.W.rowSlice(k);
            }

            const ncol = sim.wfn.W.ncol();

            for (0..ncol) |j| {
                const R_j = self.R.?.rowSlice(j);

                inline for (0..N) |k| {
                    w_col[k] = rows[k][j];
                }

                inline for (0..N) |i| {
                    var sum_re: T = 0;
                    var sum_im: T = 0;

                    inline for (0..N) |k| {
                        const r, const w = .{ R_j[i * N + k], w_col[k] };

                        sum_re += r.re * w.re - r.im * w.im;
                        sum_im += r.re * w.im + r.im * w.re;
                    }

                    rows[i][j] = Complex(T).init(sum_re, sum_im);
                }
            }
        }

        /// General potential propagation in coordinate space for arbitrary numbers of electronic states.
        fn applyRSlow(self: @This(), sim: anytype, capopt: anytype, t: T, gpa: Allocator) !void {
            const nstate = sim.wfn.W.nrow();

            var w_col = try gpa.alloc(Complex(T), nstate);
            defer gpa.free(w_col);

            for (0..sim.wfn.W.ncol()) |j| {
                const R_j = try self.getR(sim.wfn_kpgrids, sim.hams, sim.epoten, t, capopt, j);

                for (0..nstate) |k| {
                    w_col[k] = sim.wfn.W.at(k, j);
                }

                for (0..nstate) |i| {
                    var sum = Complex(T).init(0, 0);

                    for (0..nstate) |k| {
                        sum = sum.add(R_j[i * nstate + k].mul(w_col[k]));
                    }

                    sim.wfn.W.ptr(i, j).* = sum;
                }
            }
        }
    };
}
