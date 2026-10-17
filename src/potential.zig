//! Defines potential energy surfaces (PES) including ab initio, harmonic, Tully, vibronic, and custom models.

const std = @import("std");

const libint = @import("cimport.zig").libint;

const Allocator = std.mem.Allocator;

const ConfigurationInteractionOptions = @import("configuration_interaction.zig").Options;
const Expression = @import("exprtk.zig").Expression;
const HartreeFockOptions = @import("hartree_fock.zig").Options;
const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const MollerPlessetOptions = @import("moller_plesset.zig").Options;
const Value = @import("value.zig").Value;

const alignPhases = @import("configuration_interaction.zig").alignPhases;
const ci_runFromSystem = @import("configuration_interaction.zig").runFromSystem;
const exportIfBuiltin = @import("molecular_integrals.zig").exportIfBuiltin;
const generateDets = @import("configuration_interaction.zig").generateDets;
const generatePairs = @import("configuration_interaction.zig").generatePairs;
const hf_runFromSystem = @import("hartree_fock.zig").runFromSystem;
const isDual = @import("value.zig").isDual;
const mp_runFromSystem = @import("moller_plesset.zig").runFromSystem;
const readMatrix = @import("read_write.zig").readMatrix;
const stripSecondaryOptions = @import("nuclear_derivative.zig").stripSecondaryOptions;

const CM2EV = @import("constant.zig").CM2EV;
const EV2AU = @import("constant.zig").EV2AU;

/// Parameter union for potential energy surfaces supporting ab initio, harmonic, coupling, and Tully models.
pub const Options = union(enum) {
    ab_initio: AbInitioOptions,
    custom: CustomOptions,
    file: FileOptions,
    harmonic: HarmonicOptions,
    henon_heiles: HenonHeilesOptions,
    jahn_teller: JahnTellerOptions,
    lvc: LvcOptions,
    morse: MorseOptions,
    time_linear: TimeLinearOptions,
    tully_1: Tully1Options,
    tully_2: Tully2Options,
    tully_3: Tully3Options,
};

/// Ab initio electronic structure method configuration and trajectory output options.
pub const AbInitioOptions = struct {
    method: MethodOptions,
    write: AbInitioWrite = .{},
};

/// Output file destinations for ab initio potential trajectories.
pub const AbInitioWrite = struct {
    kinetic_energy: ?[]const u8 = null,
    momentum_trajectory: ?[]const u8 = null,
    population: ?[]const u8 = null,
    position_trajectory: ?[]const u8 = null,
    potential_energy: ?[]const u8 = null,
    state_coefficient: ?[]const u8 = null,
    state_potential_energy: ?[]const u8 = null,
    temperature: ?[]const u8 = null,
    total_energy: ?[]const u8 = null,

    /// Returns true if at least one trajectory property is configured for output.
    pub fn any(self: @This()) bool {
        inline for (@typeInfo(@This()).@"struct".field_names) |name| {
            if (@field(self, name) != null) return true;
        }

        return false;
    }
};

/// Available electronic structure methods for on-the-fly potential evaluation.
pub const MethodOptions = union(enum) {
    configuration_interaction: ConfigurationInteractionOptions,
    hartree_fock: HartreeFockOptions,
    moller_plesset: MollerPlessetOptions,
};

/// User-defined analytical potential matrix elements and coordinate dimension.
const CustomOptions = struct {
    matrix: []const []const []const u8,
    ndim: u32,
    time_dependent: bool = false,
};

/// Grid interpolation configuration loaded from an external potential file.
const FileOptions = struct {
    ndim: u32,
    path: []const u8,
};

/// Harmonic oscillator force constants along uncoupled Cartesian coordinates.
const HarmonicOptions = struct {
    k: []const f64 = &.{1},
};

/// Quadratic and cubic coupling parameters for the Henon-Heiles potential.
const HenonHeilesOptions = struct {
    k: f64 = 1,
    l: f64 = 0.1,
};

/// Linear and quadratic vibronic coupling parameters for E x e Jahn-Teller.
const JahnTellerOptions = struct {
    g: f64 = 1,
    k: f64 = 1,
};

/// Linear vibronic coupling Hamiltonian parameters and mode frequencies.
const LvcOptions = struct {
    excitation_energies: []const f64,
    frequencies: []const f64,
    kappa: []const []const f64,
    lambda: []const []const []const f64,
};

/// Dissociation energy, width parameter, and equilibrium distance for Morse potential.
const MorseOptions = struct {
    D: f64 = 0.2,
    a: f64 = 1,
    r0: f64 = 0,
};

/// Coupling strength and slope for time-dependent linear curve crossing.
const TimeLinearOptions = struct {
    a: f64 = 10,
    g: f64 = 2,
};

/// Model parameters for Tully 1 simple avoided crossing diabatic potential.
const Tully1Options = struct {
    A: f64 = 0.010,
    B: f64 = 1.600,
    C: f64 = 0.005,
    D: f64 = 1.000,
};

/// Model parameters for Tully 2 dual avoided crossing diabatic potential.
const Tully2Options = struct {
    A: f64 = 0.10,
    B: f64 = 0.28,
    C: f64 = 0.015,
    D: f64 = 0.06,
    E: f64 = 0.05,
};

/// Model parameters for Tully 3 extended coupling with reflection diabatic potential.
const Tully3Options = struct {
    A: f64 = 6.0e-4,
    B: f64 = 0.10,
    C: f64 = 0.90,
};

/// Returns a generic union representing a potential energy surface (PES) with coordinate and time evaluations.
pub fn Potential(comptime T: type) type {
    return union(enum) {
        ab_initio: AbInitio(T),
        custom: Custom(T),
        file: File(T),
        harmonic: Harmonic(T),
        henon_heiles: HenonHeiles(T),
        jahn_teller: JahnTeller(T),
        lvc: Lvc(T),
        morse: Morse(T),
        time_linear: TimeLinear(T),
        tully_1: Tully1(T),
        tully_2: Tully2(T),
        tully_3: Tully3(T),

        /// Initializes the selected potential energy surface based on configuration options.
        pub fn init(io: std.Io, options: Options, allocator: Allocator) !@This() {
            return switch (options) {
                .ab_initio => |f| .{ .ab_initio = try AbInitio(T).init(io, f, allocator) },
                .custom => |f| .{ .custom = try Custom(T).init(f.ndim, f.matrix, f.time_dependent, allocator) },
                .file => |f| .{ .file = try File(T).init(f.ndim, f.path, io, allocator) },
                .harmonic => |f| .{ .harmonic = Harmonic(T).init(f.k) },
                .henon_heiles => |f| .{ .henon_heiles = HenonHeiles(T).init(f.k, f.l) },
                .jahn_teller => |f| .{ .jahn_teller = JahnTeller(T).init(f.k, f.g) },
                .lvc => |f| .{ .lvc = Lvc(T).init(f.frequencies, f.excitation_energies, f.kappa, f.lambda) },
                .morse => |f| .{ .morse = Morse(T).init(f.D, f.a, f.r0) },
                .time_linear => |f| .{ .time_linear = TimeLinear(T).init(f.a, f.g) },
                .tully_1 => |f| .{ .tully_1 = Tully1(T).init(f.A, f.B, f.C, f.D) },
                .tully_2 => |f| .{ .tully_2 = Tully2(T).init(f.A, f.B, f.C, f.D, f.E) },
                .tully_3 => |f| .{ .tully_3 = Tully3(T).init(f.A, f.B, f.C) },
            };
        }

        /// Deallocates memory associated with the potential energy surface representation.
        pub fn deinit(self: *@This(), allocator: Allocator) void {
            switch (self.*) {
                inline .ab_initio, .custom, .file => |*pot| pot.deinit(allocator),

                inline else => {},
            }
        }

        /// Evaluates the potential energy matrix elements at coordinates r and time t.
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, t: U) !void {
            std.debug.assert(V.len == self.nstate() * self.nstate());

            std.debug.assert(r.len == self.ndim());

            switch (self) {
                inline else => |field| try field.eval(U, V, r, t),
            }
        }

        /// Evaluates the potential energy matrix elements for a batch of coordinate coordinates.
        pub fn evalBatch(self: @This(), comptime U: type, V: *Matrix(U), r: Matrix(U), t: T) !void {
            const t_val = Value(U).fromFloat(t);

            switch (self) {
                inline else => |field| {
                    for (0..r.nrow()) |i| try field.eval(U, V.rowSlice(i), r.rowSlice(i), t_val.val);
                },
            }
        }

        /// Returns true if the potential has explicit time dependence.
        pub fn isTd(self: @This()) bool {
            return switch (self) {
                inline else => |field| field.isTd(),
            };
        }

        /// Returns the number of spatial dimensions of the coordinates.
        pub fn ndim(self: @This()) usize {
            return switch (self) {
                inline else => |field| field.ndim(),
            };
        }

        /// Returns the number of electronic states in the potential representation.
        pub fn nstate(self: @This()) usize {
            return switch (self) {
                inline else => |field| field.nstate(),
            };
        }
    };
}

/// Returns an ab initio potential energy surface evaluating electronic energies and nuclear gradients on the fly.
fn AbInitio(comptime T: type) type {
    return struct {
        gpa: Allocator,
        inpout: std.Io,

        options: AbInitioOptions,
        msys: MolecularSystem(T),

        Pg: *std.ArrayList(?Matrix(T)),

        C_ci: *std.ArrayList(?Matrix(T)),
        C_mo: *std.ArrayList(?Matrix(T)),

        /// Initializes the molecular system and electronic structure method configuration.
        pub fn init(io: std.Io, options: AbInitioOptions, gpa: Allocator) !@This() {
            var clean_opt = options;

            clean_opt.method = stripSecondaryOptions(options.method, .{ .gradient = true, .nacv = true });

            const sys_file, const basis_name, const charge, const multi = switch (clean_opt.method) {
                .hartree_fock => |opt| .{ opt.system, opt.basis, opt.charge, opt.multiplicity },

                inline else => |opt| .{
                    opt.hartree_fock.system,
                    opt.hartree_fock.basis,
                    opt.hartree_fock.charge,
                    opt.hartree_fock.multiplicity,
                },
            };

            const has_grad = switch (clean_opt.method) {
                inline else => |opt| opt.gradient != null,
            };

            if (!has_grad) {
                std.log.err("GRADIENT CALCULATION IS MANDATORY FOR AB INITIO POTENTIAL", .{});

                return error.InvalidInput;
            }

            const basis_path = try exportIfBuiltin(io, basis_name, gpa);

            defer if (std.mem.startsWith(u8, basis_name, "builtin:")) {
                std.Io.Dir.cwd().deleteFile(io, basis_path) catch {};
            };

            var msys = try MolecularSystem(T).init(io, sys_file, basis_path, charge, multi, gpa);
            errdefer msys.deinit(gpa);

            if (std.mem.startsWith(u8, basis_name, "builtin:")) {
                try std.Io.Dir.cwd().deleteFile(io, basis_path);
            }

            const C_ci_list = try gpa.create(std.ArrayList(?Matrix(T)));
            errdefer gpa.destroy(C_ci_list);

            C_ci_list.* = .empty;

            const C_mo_list = try gpa.create(std.ArrayList(?Matrix(T)));
            errdefer gpa.destroy(C_mo_list);

            C_mo_list.* = .empty;

            const Pg_list = try gpa.create(std.ArrayList(?Matrix(T)));
            errdefer gpa.destroy(Pg_list);

            Pg_list.* = .empty;

            return .{
                .options = clean_opt,

                .inpout = io,
                .msys = msys,

                .Pg = Pg_list,

                .C_ci = C_ci_list,
                .C_mo = C_mo_list,

                .gpa = gpa,
            };
        }

        /// Deallocates molecular system resources and cached density matrix.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            for (self.C_mo.items) |*c| {
                if (c.*) |*m| m.deinit(gpa);
            }

            self.C_mo.deinit(gpa);

            gpa.destroy(self.C_mo);

            for (self.C_ci.items) |*c| {
                if (c.*) |*m| m.deinit(gpa);
            }

            self.C_ci.deinit(gpa);

            gpa.destroy(self.C_ci);

            for (self.Pg.items) |*p| {
                if (p.*) |*m| m.deinit(gpa);
            }

            self.Pg.deinit(gpa);

            gpa.destroy(self.Pg);

            self.msys.deinit(gpa);
        }

        /// Throws an error as direct potential evaluation without trajectory gradients is unsupported.
        pub fn eval(_: @This(), comptime _: type, _: anytype, _: anytype, _: anytype) !void {
            std.log.err("DIRECT EVALUATION IS NOT SUPPORTED FOR AB INITIO POTENTIAL", .{});

            return error.InvalidInput;
        }

        /// Evaluates the electronic energy and Cartesian nuclear gradients at coordinates r.
        pub fn evalGrad(self: @This(), V: []T, grad_V: []T, r: []const T, active_state: ?usize, traj: usize, eval_nacv: bool) !void {
            @memcpy(self.msys.coors, r);

            libint.libint_update_coords(self.msys.ptr, self.msys.coors.ptr);

            switch (self.options.method) {
                .hartree_fock => |opt| try self.evalGradSingle(V, grad_V, opt, hf_runFromSystem, traj),
                .moller_plesset => |opt| try self.evalGradSingle(V, grad_V, opt, mp_runFromSystem, traj),
                .configuration_interaction => |opt| try self.evalGradCi(V, grad_V, opt, active_state, traj, eval_nacv),
            }
        }

        /// Returns false as this potential is independent of time.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        /// Returns the Cartesian coordinate dimensionality of the molecular system.
        pub fn ndim(self: @This()) usize {
            return 3 * self.msys.atoms.len;
        }

        /// Returns the number of electronic states in the ab initio calculation.
        pub fn nstate(self: @This()) usize {
            return switch (self.options.method) {
                .configuration_interaction => |opt| opt.nstate,

                inline else => 1,
            };
        }

        /// Updates potential values and gradients for all batch trajectory coordinate rows.
        pub fn updateGradients(self: @This(), V: *Matrix(T), grad_V: *Matrix(T), r: Matrix(T), states: ?[]const usize) !void {
            for (0..r.nrow()) |i| {
                const s = if (states) |st| st[i] else null;

                try self.evalGrad(V.rowSlice(i), grad_V.rowSlice(i), r.rowSlice(i), s, i, true);
            }
        }

        /// Updates potential value and gradient for a specific trajectory row and active electronic state.
        pub fn updateTrajectoryGradient(self: @This(), V: *Matrix(T), grad_V: *Matrix(T), r: Matrix(T), i: usize, active_state: usize) !void {
            try self.evalGrad(V.rowSlice(i), grad_V.rowSlice(i), r.rowSlice(i), active_state, i, false);
        }

        /// Evaluates CI electronic energies, state nuclear gradients, and nonadiabatic couplings.
        fn evalGradCi(self: @This(), V: []T, grad_V: []T, opt: ConfigurationInteractionOptions, active_state: ?usize, traj: usize, eval_nacv: bool) !void {
            var run_opt = opt;

            if (!eval_nacv) {
                run_opt.nacv = null;
            }

            const n = self.nstate();

            var active_state_arr: [1]u32 = undefined;

            if (active_state) |s| {
                active_state_arr[0] = @intCast(s);
            }

            var all_states: ?[]u32 = null;
            defer if (all_states) |s| self.gpa.free(s);

            const grad_states = if (active_state != null) &active_state_arr else blk: {
                all_states = try self.gpa.alloc(u32, n);

                for (0..n) |k| all_states.?[k] = @intCast(k);

                break :blk all_states.?;
            };

            if (run_opt.gradient) |*g| switch (g.*) {
                .analytic => |*a| a.states = grad_states,
                .numeric => |*num| num.states = grad_states,
            };

            var auto_pairs: ?[][2]u32 = null;
            defer if (auto_pairs) |p| self.gpa.free(p);

            if (run_opt.nacv) |*nacvopt| switch (nacvopt.*) {
                inline else => |*num| {
                    auto_pairs = try generatePairs(n, self.gpa);

                    num.pairs = auto_pairs;
                },
            };

            while (self.Pg.items.len <= traj) try self.Pg.append(self.gpa, null);

            var res = try ci_runFromSystem(T, self.inpout, run_opt, @constCast(&self.msys), self.Pg.items[traj], false, self.gpa);
            defer res.deinit(self.gpa);

            if (self.Pg.items[traj] == null) {
                self.Pg.items[traj] = try Matrix(T).init(res.hartree_fock.P.nrow(), res.hartree_fock.P.ncol(), self.gpa);
            }

            @memcpy(self.Pg.items[traj].?.data, res.hartree_fock.P.data);

            while (self.C_mo.items.len <= traj) try self.C_mo.append(self.gpa, null);
            while (self.C_ci.items.len <= traj) try self.C_ci.append(self.gpa, null);

            const flipped = try self.gpa.alloc(bool, res.C.ncol());
            defer self.gpa.free(flipped);

            @memset(flipped, false);

            const gen = opt.hartree_fock.generalized;

            var dets = try generateDets(self.msys.nel, 2 * self.msys.nbf, opt.excitations, self.gpa);

            defer {
                for (0..dets.items.len) |idx| self.gpa.free(dets.items[idx]);

                dets.deinit(self.gpa);
            }

            const hf_C, const ci_C = .{ &res.hartree_fock.C, &res.C };

            alignPhases(T, hf_C, ci_C, self.C_mo.items[traj], self.C_ci.items[traj], dets.items, gen, flipped);

            if (self.C_mo.items[traj] == null) {
                self.C_mo.items[traj] = try Matrix(T).init(hf_C.nrow(), hf_C.ncol(), self.gpa);
            }

            @memcpy(self.C_mo.items[traj].?.data, res.hartree_fock.C.data);

            if (self.C_ci.items[traj] == null) {
                self.C_ci.items[traj] = try Matrix(T).init(res.C.nrow(), res.C.ncol(), self.gpa);
            }

            @memcpy(self.C_ci.items[traj].?.data, res.C.data);

            @memset(V[0 .. n * n], 0);

            if (eval_nacv) {
                @memset(grad_V[0 .. self.ndim() * n * n], 0);
            }

            for (0..n) |k| {
                V[k * n + k] = res.energy[k];
            }

            if (res.grad.len > 0) {
                if (active_state) |s| {
                    for (0..self.ndim()) |j| {
                        grad_V[j * n * n + s * n + s] = res.grad[0].data[j];
                    }
                }

                if (active_state == null) {
                    for (0..n) |k| for (0..self.ndim()) |j| {
                        grad_V[j * n * n + k * n + k] = res.grad[k].data[j];
                    };
                }
            }

            if (auto_pairs) |pairs| {
                for (pairs, 0..) |pair, idx| {
                    const p0, const p1 = .{ pair[0], pair[1] };

                    const sign: T = if (flipped[p0] != flipped[p1]) -1 else 1;

                    for (0..self.ndim()) |j| {
                        const val = sign * res.nacv[idx].data[j];

                        grad_V[j * n * n + p0 * n + p1] = val;
                        grad_V[j * n * n + p1 * n + p0] = -val;
                    }
                }
            }
        }

        /// Evaluates electronic energy and nuclear gradients for single-reference electronic structure methods.
        fn evalGradSingle(self: @This(), V: []T, grad_V: []T, opt: anytype, comptime runFn: anytype, traj: usize) !void {
            while (self.Pg.items.len <= traj) try self.Pg.append(self.gpa, null);

            var res = try runFn(T, self.inpout, opt, @constCast(&self.msys), self.Pg.items[traj], false, self.gpa);
            defer res.deinit(self.gpa);

            const P = if (@hasField(@TypeOf(res), "hartree_fock")) res.hartree_fock.P else res.P;

            if (self.Pg.items[traj] == null) {
                self.Pg.items[traj] = try Matrix(T).init(P.nrow(), P.ncol(), self.gpa);
            }

            @memcpy(self.Pg.items[traj].?.data, P.data);

            V[0] = res.energy[0];

            @memcpy(grad_V[0..self.ndim()], res.grad[0].data[0..self.ndim()]);
        }
    };
}

/// Returns a custom potential type parsed from analytical string expressions.
fn Custom(comptime T: type) type {
    return struct {
        expressions: []const Expression(T),

        is_td: bool,
        r_vals: []T,

        /// Parses custom potential matrix expressions for a given coordinate dimension.
        pub fn init(dim: usize, matrix: []const []const []const u8, is_td: bool, allocator: Allocator) !@This() {
            for (0..matrix.len) |i| {
                if (matrix[i].len != matrix.len) return error.NonSquareMatrix;
            }

            const exprs = try allocator.alloc(Expression(T), matrix.len * matrix.len);
            errdefer allocator.free(exprs);

            for (0..matrix.len) |i| for (0..matrix.len) |j| {
                exprs[i * matrix.len + j] = Expression(T).init(matrix[i][j], dim) catch |err| {
                    for (0..(i * matrix.len + j)) |k| {
                        exprs[k].deinit();
                    }

                    return err;
                };
            };

            errdefer for (0..matrix.len * matrix.len) |i| {
                exprs[i].deinit();
            };

            return .{ .expressions = exprs, .r_vals = try allocator.alloc(T, dim), .is_td = is_td };
        }

        /// Frees parsed expression resources and internal coordinate buffer.
        pub fn deinit(self: *@This(), allocator: Allocator) void {
            allocator.free(self.r_vals);

            for (self.expressions) |expr| {
                expr.deinit();
            }

            allocator.free(self.expressions);
        }

        /// Evaluates custom potentials by parsing expressions using ExprTk.
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, t: U) !void {
            if (comptime U == T) for (0..self.nstate()) |i| for (0..self.nstate()) |j| {
                V[i * self.nstate() + j] = self.expressions[i * self.nstate() + j].evaluate_d0(r, t);
            };

            if (comptime isDual(U)) {
                for (0..r.len) |i| {
                    self.r_vals[i] = r[i].val;
                }

                for (0..self.nstate()) |i| for (0..self.nstate()) |j| {
                    const expr, var der: T = .{ self.expressions[i * self.nstate() + j], 0 };

                    for (0..r.len) |k| {
                        der += expr.evaluate_d1(self.r_vals, t.val, k) * r[k].der;
                    }

                    if (t.der != 0) {
                        der += expr.evaluate_d1(self.r_vals, t.val, r.len) * t.der;
                    }

                    V[i * self.nstate() + j] = U.init(expr.evaluate_d0(self.r_vals, t.val), der);
                };
            }

            if (comptime U != T and !isDual(U)) {
                @compileError("UNSUPPORTED NUMBER TYPE FOR CUSTOM POTENTIAL");
            }
        }

        /// Returns true if the custom potential expression explicitly depends on time.
        pub fn isTd(self: @This()) bool {
            return self.is_td;
        }

        /// Returns the spatial coordinate dimensionality of the custom potential.
        pub fn ndim(self: @This()) usize {
            return self.r_vals.len;
        }

        /// Returns the number of coupled electronic states in the custom potential matrix.
        pub fn nstate(self: @This()) usize {
            return std.math.sqrt(self.expressions.len);
        }
    };
}

/// Returns a potential type interpolated from grid data stored in a file.
fn File(comptime T: type) type {
    return struct {
        U: Matrix(T),
        ndims: usize,

        /// Reads the potential grid data from the given file path.
        pub fn init(ndims: usize, path: []const u8, io: std.Io, gpa: Allocator) !@This() {
            return .{ .U = try readMatrix(T, io, path, gpa), .ndims = ndims };
        }

        /// Deallocates the stored potential grid matrix.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.U.deinit(gpa);
        }

        /// Evaluates the potential via multilinear interpolation of grid values.
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) !void {
            for (0..self.nstate()) |i| for (0..self.nstate()) |j| {
                const col = i * self.nstate() + j;

                V[col] = lerp(T, U, self.U, self.ndim() + col, r);
            };
        }

        /// Returns false as the grid file potential is independent of time.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        /// Returns the spatial coordinate dimensionality of the grid.
        pub fn ndim(self: @This()) usize {
            return self.ndims;
        }

        /// Returns the number of coupled electronic states in the interpolated potential.
        pub fn nstate(self: @This()) usize {
            return std.math.sqrt(self.U.ncol() - self.ndim());
        }
    };
}

/// Returns a harmonic potential energy surface type parameterized by force constants k.
fn Harmonic(comptime T: type) type {
    return struct {
        k: []const T,

        /// Initializes a harmonic potential with specified force constants.
        pub fn init(k: []const T) @This() {
            return .{ .k = k };
        }

        /// Evaluates the harmonic potential energy: V = 0.5 * sum( k_i * r_i^2 ).
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) !void {
            var sum = Value(U).fromFloat(0);

            for (0..r.len) |i| {
                sum = sum.add(Value(U).init(r[i]).mul(Value(U).init(r[i])).muls(0.5 * self.k[i]));
            }

            V[0] = sum.val;
        }

        /// Returns false as the harmonic potential is independent of time.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        /// Returns the number of harmonic oscillator dimensions.
        pub fn ndim(self: @This()) usize {
            return self.k.len;
        }

        /// Returns the single ground state of the uncoupled harmonic potential.
        pub fn nstate(_: @This()) usize {
            return 1;
        }
    };
}

/// Returns a Henon-Heiles potential energy surface type parameterized by frequency and anharmonicity.
fn HenonHeiles(comptime T: type) type {
    return struct {
        k: T,
        l: T,

        /// Initializes a Henon-Heiles potential with frequency and anharmonicity parameters.
        pub fn init(k: T, l: T) @This() {
            return .{ .k = k, .l = l };
        }

        /// Evaluates the Henon-Heiles potential energy: V = 0.5 * omg^2 * (x^2 + y^2) + lmb * (x^2 * y - y^3 / 3).
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) !void {
            std.debug.assert(r.len == 2);

            const x = Value(U).init(r[0]);
            const y = Value(U).init(r[1]);

            const k = Value(U).fromFloat(self.k);
            const l = Value(U).fromFloat(self.l);

            const V0 = k.muls(0.5).mul(x.mul(x).add(y.mul(y)));
            const V1 = l.mul(x.mul(x).mul(y).sub(y.mul(y).mul(y).divs(3)));

            V[0] = V0.add(V1).val;
        }

        /// Returns false as this potential is independent of time.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        /// Returns the two spatial coordinate dimensions of the Henon-Heiles system.
        pub fn ndim(_: @This()) usize {
            return 2;
        }

        /// Returns the single electronic state of the Henon-Heiles potential.
        pub fn nstate(_: @This()) usize {
            return 1;
        }
    };
}

/// Returns a Jahn-Teller potential energy surface type parameterized by frequency and coupling strength.
fn JahnTeller(comptime T: type) type {
    return struct {
        k: T,
        g: T,

        /// Initializes a Jahn-Teller potential with frequency and coupling strength parameters.
        pub fn init(k: T, g: T) @This() {
            return .{ .k = k, .g = g };
        }

        /// Evaluates the Jahn-Teller potential energy matrix: V = 0.5 * k * (x^2 + y^2) * I + g * [[x, y], [y, -x]].
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) !void {
            std.debug.assert(r.len == 2);

            const x = Value(U).init(r[0]);
            const y = Value(U).init(r[1]);

            const k = Value(U).fromFloat(self.k);
            const g = Value(U).fromFloat(self.g);

            const V00 = k.muls(0.5).mul(x.mul(x).add(y.mul(y))).add(g.mul(x));
            const V01 = g.mul(y);
            const V11 = k.muls(0.5).mul(x.mul(x).add(y.mul(y))).sub(g.mul(x));

            V[0] = V00.val;
            V[1] = V01.val;
            V[2] = V01.val;
            V[3] = V11.val;
        }

        /// Returns false as the Jahn-Teller potential is independent of time.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        /// Returns the two degenerate normal mode coordinate dimensions.
        pub fn ndim(_: @This()) usize {
            return 2;
        }

        /// Returns the two electronic states coupled by the Jahn-Teller interaction.
        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}

/// Returns a linear vibronic coupling (LVC) potential energy surface type parameterized in eV and cm^-1.
fn Lvc(comptime T: type) type {
    return struct {
        frequencies: []const T,
        kap: []const []const T,

        excitation_energies: []const T,
        lmb: []const []const []const T,

        /// Initializes the linear vibronic coupling potential with frequencies, excitations, and coupling constants.
        pub fn init(frequencies: []const T, exc_en: []const T, kappa: []const []const T, lambda: []const []const []const T) @This() {
            std.debug.assert(kappa.len == exc_en.len);

            for (kappa) |row| {
                std.debug.assert(row.len == frequencies.len);
            }

            std.debug.assert(lambda.len == exc_en.len);

            for (lambda) |row| {
                std.debug.assert(row.len == exc_en.len);

                for (row) |col| {
                    std.debug.assert(col.len == frequencies.len);
                }
            }

            return .{ .frequencies = frequencies, .excitation_energies = exc_en, .kap = kappa, .lmb = lambda };
        }

        /// Evaluates the vibronic coupling potential matrix in Hartree at dimensionless normal coordinates r.
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) !void {
            std.debug.assert(V.len == self.nstate() * self.nstate());

            var v0 = Value(U).fromFloat(0);

            for (0..self.frequencies.len) |i| {
                const q = Value(U).init(r[i]);

                v0 = v0.add(q.mul(q).muls(0.5 * self.frequencies[i] * CM2EV));
            }

            for (0..self.excitation_energies.len) |n| {
                var term = Value(U).fromFloat(self.excitation_energies[n]);

                for (0..self.frequencies.len) |i| {
                    term = term.add(Value(U).init(r[i]).muls(self.kap[n][i]));
                }

                V[n * self.nstate() + n] = v0.add(term).val;
            }

            for (0..self.excitation_energies.len) |m| for (m + 1..self.excitation_energies.len) |n| {
                var term = Value(U).fromFloat(0);

                for (0..self.frequencies.len) |i| {
                    term = term.add(Value(U).init(r[i]).muls(self.lmb[m][n][i]));
                }

                V[m * self.nstate() + n] = term.val;
                V[n * self.nstate() + m] = term.val;
            };

            for (0..V.len) |i| {
                V[i] = Value(U).init(V[i]).muls(EV2AU).val;
            }
        }

        /// Returns false as this vibronic potential has no explicit time dependence.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        /// Returns the number of normal modes (nuclear degrees of freedom) in the potential.
        pub fn ndim(self: @This()) usize {
            return self.frequencies.len;
        }

        /// Returns the number of electronic states in the non-adiabatic potential representation.
        pub fn nstate(self: @This()) usize {
            return self.excitation_energies.len;
        }
    };
}

/// Returns a Morse potential energy surface type describing anharmonic bond stretching.
fn Morse(comptime T: type) type {
    return struct {
        D: T,
        a: T,

        r0: T,

        /// Initializes a Morse potential with dissociation energy D, stiffness a, and equilibrium distance r0.
        pub fn init(D: T, a: T, r0: T) @This() {
            return .{ .D = D, .a = a, .r0 = r0 };
        }

        /// Evaluates the Morse potential energy: V = D * (1 - exp(-a * (r - r0)))^2.
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) !void {
            const r0 = Value(U).init(r[0]);

            const D = Value(U).fromFloat(self.D);
            const a = Value(U).fromFloat(self.a);

            const r_eq = Value(U).fromFloat(self.r0);

            const factor = Value(U).fromFloat(1).sub(r0.sub(r_eq).mul(a).neg().exp());

            V[0] = D.mul(factor).mul(factor).val;
        }

        /// Returns false as this potential is independent of time.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        /// Returns the single nuclear coordinate dimension of the Morse potential.
        pub fn ndim(_: @This()) usize {
            return 1;
        }

        /// Returns the single electronic ground state of the Morse potential.
        pub fn nstate(_: @This()) usize {
            return 1;
        }
    };
}

/// Returns a time-dependent linear two-state coupling potential.
fn TimeLinear(comptime T: type) type {
    return struct {
        a: T,
        g: T,

        /// Initializes the linear coupling potential with slope and coupling parameters.
        pub fn init(a: T, g: T) @This() {
            return .{ .a = a, .g = g };
        }

        /// Evaluates the two-state time-dependent linear coupling potential matrix.
        pub fn eval(self: @This(), comptime U: type, V: []U, _: []const U, t: U) !void {
            const a = Value(U).fromFloat(self.a);
            const g = Value(U).fromFloat(self.g);

            const V00 = Value(U).init(t).sub(a).mul(a);
            const V01 = g;
            const V11 = V00.neg();

            V[0] = V00.val;
            V[1] = V01.val;
            V[2] = V01.val;
            V[3] = V11.val;
        }

        /// Returns true as the linear coupling potential explicitly depends on time.
        pub fn isTd(_: @This()) bool {
            return true;
        }

        /// Returns the single nuclear coordinate dimension of the time-linear model.
        pub fn ndim(_: @This()) usize {
            return 1;
        }

        /// Returns the two coupled electronic states of the time-linear model.
        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}

/// Returns a Tully 1 simple avoided crossing potential model type.
fn Tully1(comptime T: type) type {
    return struct {
        A: T,
        B: T,
        C: T,
        D: T,

        /// Initializes Tully 1 model parameters representing simple avoided crossing.
        pub fn init(A: T, B: T, C: T, D: T) @This() {
            return .{ .A = A, .B = B, .C = C, .D = D };
        }

        /// Evaluates the Tully simple avoided crossing two-state potential matrix.
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) !void {
            const r0 = Value(U).init(r[0]);

            const A = Value(U).fromFloat(self.A);
            const B = Value(U).fromFloat(self.B);
            const C = Value(U).fromFloat(self.C);
            const D = Value(U).fromFloat(self.D);

            const V00 = r0.abs().mul(B).neg().exp().neg().adds(1).mul(A).mul(r0.sign());
            const V01 = r0.mul(r0).mul(D).neg().exp().mul(C);
            const V11 = V00.neg();

            V[0] = V00.val;
            V[1] = V01.val;
            V[2] = V01.val;
            V[3] = V11.val;
        }

        /// Returns false as this potential is independent of time.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        /// Returns the single nuclear coordinate dimension of the Tully 1 model.
        pub fn ndim(_: @This()) usize {
            return 1;
        }

        /// Returns the two coupled electronic states in the Tully 1 model.
        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}

/// Returns a Tully 2 dual avoided crossing potential model type.
fn Tully2(comptime T: type) type {
    return struct {
        A: T,
        B: T,
        C: T,
        D: T,
        E: T,

        /// Initializes Tully 2 model parameters representing dual avoided crossing.
        pub fn init(A: T, B: T, C: T, D: T, E: T) @This() {
            return .{ .A = A, .B = B, .C = C, .D = D, .E = E };
        }

        /// Evaluates the Tully dual avoided crossing two-state potential matrix.
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) !void {
            const r0 = Value(U).init(r[0]);

            const A = Value(U).fromFloat(self.A);
            const B = Value(U).fromFloat(self.B);
            const C = Value(U).fromFloat(self.C);
            const D = Value(U).fromFloat(self.D);

            const E = Value(U).fromFloat(self.E);

            const V00 = Value(U).fromFloat(0);
            const V01 = r0.mul(r0).mul(D).neg().exp().mul(C);
            const V11 = r0.mul(r0).mul(B).neg().exp().mul(A).neg().add(E);

            V[0] = V00.val;
            V[1] = V01.val;
            V[2] = V01.val;
            V[3] = V11.val;
        }

        /// Returns false as this potential is independent of time.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        /// Returns the single nuclear coordinate dimension of the model.
        pub fn ndim(_: @This()) usize {
            return 1;
        }

        /// Returns the number of electronic states in this model system.
        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}

/// Returns a Tully 3 extended coupling with reflection potential model type.
fn Tully3(comptime T: type) type {
    return struct {
        A: T,
        B: T,
        C: T,

        /// Initializes Tully 3 model parameters representing extended coupling.
        pub fn init(A: T, B: T, C: T) @This() {
            return .{ .A = A, .B = B, .C = C };
        }

        /// Evaluates the Tully extended coupling two-state potential matrix.
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) !void {
            const r0 = Value(U).init(r[0]);

            const A = Value(U).fromFloat(self.A);
            const B = Value(U).fromFloat(self.B);
            const C = Value(U).fromFloat(self.C);

            const V00 = A;
            const V01 = Value(U).fromFloat(1).sub(r0.abs().mul(C).neg().exp()).mul(r0.sign()).adds(1).mul(B);
            const V11 = A.neg();

            V[0] = V00.val;
            V[1] = V01.val;
            V[2] = V01.val;
            V[3] = V11.val;
        }

        /// Returns false as this potential is independent of time.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        /// Returns the single nuclear coordinate dimension of the model.
        pub fn ndim(_: @This()) usize {
            return 1;
        }

        /// Returns the number of electronic states in this model system.
        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}

/// Performs multilinear interpolation on a rectangular coordinate grid.
fn lerp(comptime T: type, comptime U: type, grid: Matrix(T), column: usize, r: []const U) U {
    std.debug.assert(column < grid.ncol());

    const size: usize = @round(std.math.pow(T, @as(T, @floatFromInt(grid.nrow())), 1 / @as(T, @floatFromInt(r.len))));

    var result = Value(U).fromFloat(0);

    for (0..@as(usize, 1) << @as(u5, @intCast(r.len))) |i| {
        var w, var j: usize = .{ Value(U).fromFloat(1), 0 };

        for (0..r.len) |k| {
            const stride = std.math.pow(usize, size, r.len - k - 1);

            var low: usize, var high: usize, var mid: usize = .{ 0, size, size / 2 };

            while (low < high) : (mid = (low + high) / 2) {
                const r_k = if (comptime isDual(U)) r[k].val else r[k];

                if (grid.at(mid * stride, k) <= r_k) low = mid + 1 else high = mid;
            }

            const idx = @min(@max(low, 1), size - 1);

            const x0 = grid.at((idx - 1) * stride, k);
            const x1 = grid.at((idx - 0) * stride, k);

            const wgh = Value(U).init(r[k]).sub(Value(U).fromFloat(x0)).muls(1 / (x1 - x0));

            const use_upper = ((i >> @as(u5, @intCast(r.len - k - 1))) & 1) == 1;

            const wmul = if (use_upper) wgh else Value(U).fromFloat(1).sub(wgh);

            w, j = .{ w.mul(wmul), j + (if (use_upper) idx else idx - 1) * stride };
        }

        result = result.add(w.muls(grid.at(j, column)));
    }

    return result.val;
}
