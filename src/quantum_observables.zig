//! Evaluates quantum mechanical expectation values and autocorrelation functions on spatial grids.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Grid = @import("wavepacket.zig").Grid;
const Hamiltonian = @import("wavepacket.zig").Hamiltonian;
const Potential = @import("potential.zig").Potential;
const Vector = @import("tensor.zig").Vector;
const Wavefunction = @import("wavepacket.zig").Wavefunction;

/// Computes instantaneous thermodynamic and nonadiabatic quantum expectation values.
pub fn Observables(comptime T: type) type {
    return struct {
        pos: ?Vector(T) = null,
        mom: ?Vector(T) = null,
        pop: ?Vector(T) = null,

        acf: ?Complex(T) = null,

        epot: ?T = null,
        ekin: ?T = null,
        norm: ?T = null,

        /// Computes quantum expectation values and state populations from the current wavepacket.
        pub fn init(sim: anytype, wfn0: ?Wavefunction(T), write: anytype, adia: bool, log: bool, t: T, gpa: Allocator) !@This() {
            var obs = @This(){};
            errdefer obs.deinit(gpa);

            const calc_ekin, const calc_epot = .{ write.kinetic_energy != null, write.potential_energy != null };

            const calc_norm = write.norm != null;

            var calc = .{
                .pos = log or write.position != null,
                .mom = log or write.momentum != null,
                .pop = log or write.population != null,
                .norm = log or calc_norm,
                .ekin = log or calc_ekin,
                .epot = log or calc_epot,
            };

            if (wfn0 != null) {
                obs.acf = wfn0.?.overlap(sim.wfn, sim.wfn_kpgrids);
            }

            calc.ekin = calc.ekin or write.total_energy != null;
            calc.epot = calc.epot or write.total_energy != null;

            if (calc.pos) {
                obs.pos = try sim.wfn.pos(sim.wfn_kpgrids, gpa);
            }

            if (calc.norm) {
                obs.norm = sim.wfn.norm(sim.wfn_kpgrids);
            }

            const langer = computeLangerCorrection(T, sim.wfn, sim.wfn_kpgrids, sim.hams, calc.epot or calc.ekin);

            if (calc.epot) {
                obs.epot = try sim.wfn.epot(sim.hams, sim.wfn_kpgrids, sim.epoten, t, langer);
            }

            if (calc.pop) {
                if (adia == true) {
                    obs.pop = try sim.wfn.popAdia(sim.hams, sim.wfn_kpgrids, sim.epoten, t, gpa);
                }

                if (adia != true) {
                    obs.pop = try sim.wfn.pop(sim.wfn_kpgrids, gpa);
                }

                if (obs.pop) |*pop| for (0..pop.length()) |j| {
                    pop.ptr(j).* += sim.pop_apabs.at(j);
                };
            }

            const needs_fft = calc.mom or calc.ekin;

            if (needs_fft) {
                sim.wfn.fft(-1);

                defer {
                    sim.wfn.fft(1);
                }

                if (calc.mom) {
                    obs.mom = try sim.wfn.mom(sim.wfn_kpgrids, gpa);
                }

                if (calc.ekin) {
                    obs.ekin = sim.wfn.ekin(sim.hams, sim.wfn_kpgrids, langer);
                }
            }

            return obs;
        }

        /// Deallocates arrays stored in the quantum observables struct.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.pos) |*pos| pos.deinit(gpa);
            if (self.mom) |*mom| mom.deinit(gpa);
            if (self.pop) |*pop| pop.deinit(gpa);
        }
    };
}

/// Evaluates the Langer quantum centrifugal correction for cylindrical coordinate grids.
fn computeLangerCorrection(comptime T: type, wfn: anytype, grid: anytype, hams: anytype, needed: bool) T {
    if (!hams.cylindric or !needed) return 0;

    var langer: T = 0;

    const radial_idx = grid.ncol() - 1;

    const m = hams.mass[radial_idx];

    for (0..wfn.W.nrow()) |s| for (0..wfn.W.rowSlice(s).len) |j| {
        const r = grid.getR(j, radial_idx);

        if (r != 0) {
            langer += wfn.W.at(s, j).squaredMagnitude() / (8 * m * r * r);
        }
    };

    return langer * grid.dr;
}
