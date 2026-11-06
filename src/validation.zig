//! Validates configuration parameters and physical consistency constraints across molecular simulation targets.

const std = @import("std");

const casscf = @import("casscf.zig");
const classical_dynamics = @import("classical_dynamics.zig");
const configuration_interaction = @import("configuration_interaction.zig");
const hartree_fock = @import("hartree_fock.zig");
const molecular_integrals = @import("molecular_integrals.zig");
const moller_plesset = @import("moller_plesset.zig");
const quantum_dynamics = @import("quantum_dynamics.zig");

/// Dispatches configuration validation to target-specific physical constraint checkers based on option types.
pub fn checkInvalidInput(opt: anytype) !void {
    const Opt = @TypeOf(opt);

    if (Opt == casscf.Options) {
        return checkCasscf(opt);
    }

    if (Opt == classical_dynamics.Options) {
        return checkClassicalDynamics(opt);
    }

    if (Opt == configuration_interaction.Options) {
        return checkConfigurationInteraction(opt);
    }

    if (Opt == hartree_fock.Options) {
        return checkHartreeFock(opt);
    }

    if (Opt == molecular_integrals.Options) {
        return checkMolecularIntegrals(opt);
    }

    if (Opt == moller_plesset.Options) {
        return checkMollerPlesset(opt);
    }

    if (Opt == quantum_dynamics.Options) {
        return checkQuantumDynamics(opt);
    }
}

/// Validates complete active space orbital ranges, electron counts, and multi-state coupling constraints.
fn checkCasscf(opt: casscf.Options) !void {
    if (opt.hartree_fock.generalized) {
        std.log.err("GENERALIZED HARTREE-FOCK IS NOT SUPPORTED FOR CASSCF", .{});

        return error.InvalidInput;
    }

    if (opt.hartree_fock.integral_direct) {
        std.log.err("INTEGRAL DIRECT CALCULATION IS NOT SUPPORTED FOR CASSCF", .{});

        return error.InvalidInput;
    }

    if (opt.active_electrons == 0 or opt.active_orbitals == 0) {
        std.log.err("CASSCF ACTIVE ELECTRONS AND ACTIVE ORBITALS MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.multiplicity == 0) {
        std.log.err("CASSCF MULTIPLICITY MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.spin_penalty) |pen| if (pen < 0 or std.math.isNan(pen)) {
        std.log.err("CASSCF SPIN PENALTY MUST BE NON-NEGATIVE", .{});

        return error.InvalidInput;
    };

    const nt, const ni, const ns = .{ opt.nthreads == 0, opt.iterations == 0, opt.nstate == 0 };

    if (nt or ni or opt.threshold <= 0 or opt.orbital_threshold <= 0 or ns) {
        std.log.err("INVALID THREADS, ITERATIONS, THRESHOLD, OR NSTATE", .{});

        return error.InvalidInput;
    }

    if (opt.weights) |weights| {
        if (weights.len != opt.nstate) {
            std.log.err("CASSCF WEIGHTS COUNT MUST MATCH NUMBER OF STATES", .{});

            return error.InvalidInput;
        }

        var sum: f64 = 0;

        for (weights) |w| {
            if (w < 0) {
                std.log.err("CASSCF STATE WEIGHTS MUST BE NON-NEGATIVE", .{});

                return error.InvalidInput;
            }

            sum += w;
        }

        if (sum <= 0) {
            std.log.err("CASSCF STATE WEIGHTS SUM MUST BE GREATER THAN ZERO", .{});

            return error.InvalidInput;
        }
    }

    try checkWriteCalculated(opt);

    try checkMultiState("CASSCF", opt);
}

/// Validates time steps, initial phase space distributions, masses, and potential coupling for trajectories.
fn checkClassicalDynamics(opt: classical_dynamics.Options) !void {
    if (opt.time_step <= 0) {
        std.log.err("TIME STEP MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.trajectories == 0) {
        std.log.err("NUMBER OF TRAJECTORIES MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.log_interval == 0) {
        std.log.err("LOG INTERVAL MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.mass) |mass_vec| {
        if (opt.potential == .ab_initio) {
            std.log.err("MASS VECTOR MUST NOT BE PROVIDED FOR AB INITIO POTENTIAL", .{});

            return error.InvalidInput;
        }

        for (mass_vec) |m| if (m <= 0) {
            std.log.err("MASS MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        };

        if (opt.initial_conditions.position) |pos| if (mass_vec.len != pos.len) {
            std.log.err("MASS VECTOR MUST HAVE THE SAME LENGTH AS POSITION VECTOR", .{});

            return error.InvalidInput;
        };
    }

    if (opt.mass == null and opt.potential != .ab_initio) {
        std.log.err("MASS VECTOR MUST BE PROVIDED FOR MODEL POTENTIALS", .{});

        return error.InvalidInput;
    }

    if (opt.initial_conditions.position) |pos| {
        if (pos.len == 0) {
            std.log.err("INITIAL POSITION VECTOR MUST NOT BE EMPTY", .{});

            return error.InvalidInput;
        }

        if (opt.initial_conditions.momentum) |mom| if (mom.len != pos.len) {
            std.log.err("INITIAL MOMENTUM AND POSITION VECTORS MUST HAVE THE SAME LENGTH", .{});

            return error.InvalidInput;
        };
    }

    if (opt.initial_conditions.distribution) |dist| switch (dist) {
        .maxwell_boltzmann => |mb| {
            if (std.math.isNan(mb.temperature) or mb.temperature < 0) {
                std.log.err("MAXWELL-BOLTZMANN TEMPERATURE MUST BE NON-NEGATIVE", .{});

                return error.InvalidInput;
            }
        },

        .wigner => |wig| {
            if (wig.temperature) |t| if (std.math.isNan(t) or t < 0) {
                std.log.err("WIGNER TEMPERATURE MUST BE NON-NEGATIVE", .{});

                return error.InvalidInput;
            };

            if (wig.gamma != null and wig.method != null) {
                std.log.err("GAMMA AND METHOD MUST NOT BE SPECIFIED SIMULTANEOUSLY FOR WIGNER SAMPLING", .{});

                return error.InvalidInput;
            }

            if (wig.gamma == null and wig.method == null and opt.potential != .ab_initio) {
                std.log.err("EITHER GAMMA OR METHOD MUST BE SPECIFIED FOR WIGNER SAMPLING", .{});

                return error.InvalidInput;
            }

            if (wig.gamma) |gam| {
                if (wig.temperature != null) {
                    std.log.err("TEMPERATURE IS NOT SUPPORTED FOR WIGNER SAMPLING WITH GAMMA", .{});

                    return error.InvalidInput;
                }

                if (opt.potential == .ab_initio) {
                    std.log.err("GAMMA VECTOR IS NOT SUPPORTED FOR AB INITIO POTENTIAL", .{});

                    return error.InvalidInput;
                }

                if (gam.len == 0) {
                    std.log.err("WIGNER GAMMA VECTOR MUST NOT BE EMPTY", .{});

                    return error.InvalidInput;
                }

                for (gam) |g| if (std.math.isNan(g) or g <= 0) {
                    std.log.err("WIGNER GAMMA VALUES MUST BE GREATER THAN 0", .{});

                    return error.InvalidInput;
                };

                if (opt.initial_conditions.position) |pos| if (gam.len != pos.len) {
                    std.log.err("WIGNER GAMMA VECTOR MUST HAVE THE SAME LENGTH AS POSITION VECTOR", .{});

                    return error.InvalidInput;
                };
            }

            if (wig.method != null or opt.potential == .ab_initio) {
                if (opt.initial_conditions.position != null) {
                    std.log.err("INITIAL POSITION MUST NOT BE SPECIFIED FOR AB INITIO WIGNER SAMPLING", .{});

                    return error.InvalidInput;
                }

                if (opt.initial_conditions.momentum != null) {
                    std.log.err("INITIAL MOMENTUM MUST NOT BE SPECIFIED FOR AB INITIO WIGNER SAMPLING", .{});

                    return error.InvalidInput;
                }
            }

            if (wig.method) |mopt| {
                if (opt.potential != .ab_initio) {
                    std.log.err("AB INITIO METHOD FOR WIGNER SAMPLING IS NOT SUPPORTED FOR MODEL POTENTIALS", .{});

                    return error.InvalidInput;
                }

                const has_hessian = switch (mopt) {
                    inline else => |m| m.hessian != null,
                };

                if (!has_hessian) {
                    std.log.err("HESSIAN MUST BE CONFIGURED IN ELECTRONIC STRUCTURE METHOD FOR WIGNER SAMPLING", .{});

                    return error.InvalidInput;
                }

                const has_gradient = switch (mopt) {
                    inline else => |m| m.gradient != null,
                };

                if (has_gradient) {
                    std.log.err("GRADIENT MUST NOT BE SPECIFIED IN ELECTRONIC STRUCTURE FOR WIGNER SAMPLING", .{});

                    return error.InvalidInput;
                }

                if (mopt == .configuration_interaction) {
                    const ci = mopt.configuration_interaction;

                    if (ci.nacv != null) {
                        std.log.err("NACV MUST NOT BE SPECIFIED IN CONFIGURATION INTERACTION FOR WIGNER SAMPLING", .{});

                        return error.InvalidInput;
                    }

                    if (ci.hartree_fock.gradient != null) {
                        std.log.err("GRADIENT MUST NOT BE SPECIFIED IN HARTREE-FOCK FOR WIGNER SAMPLING", .{});

                        return error.InvalidInput;
                    }

                    if (ci.hessian) |hess| {
                        const hess_states = switch (hess) {
                            inline else => |h| h.states,
                        };

                        if (hess_states != null) {
                            std.log.err("HESSIAN STATE MUST NOT BE SPECIFIED FOR WIGNER SAMPLING", .{});

                            return error.InvalidInput;
                        }
                    }
                }
            }
        },
    };

    if (opt.initial_conditions.position == null and opt.potential != .ab_initio) {
        std.log.err("INITIAL POSITION VECTOR MUST BE PROVIDED FOR MODEL POTENTIALS", .{});

        return error.InvalidInput;
    }

    if (opt.potential == .ab_initio) {
        switch (opt.potential.ab_initio) {
            inline else => |m| {
                if (m.optimize != null) {
                    std.log.err("GEOMETRY OPTIMIZATION IS NOT SUPPORTED FOR AB INITIO POTENTIAL", .{});

                    return error.InvalidInput;
                }

                if (m.frequency != null) {
                    std.log.err("FREQUENCY CALCULATION IS NOT SUPPORTED FOR AB INITIO POTENTIAL", .{});

                    return error.InvalidInput;
                }

                if (m.hessian != null) {
                    std.log.err("HESSIAN CALCULATION IS NOT SUPPORTED FOR AB INITIO POTENTIAL", .{});

                    return error.InvalidInput;
                }
            },
        }

        const has_grad = switch (opt.potential.ab_initio) {
            inline else => |m| m.gradient != null,
        };

        if (!has_grad) {
            std.log.err("GRADIENT CALCULATION IS MANDATORY FOR AB INITIO POTENTIAL", .{});

            return error.InvalidInput;
        }
    }

    if (opt.potential == .ab_initio and !opt.adiabatic) {
        std.log.err("DIABATIC SIMULATION IS NOT SUPPORTED FOR AB INITIO POTENTIAL", .{});

        return error.InvalidInput;
    }

    if (opt.potential == .ab_initio and opt.nonadiabatic != null) {
        if (opt.potential.ab_initio == .configuration_interaction) {
            if (opt.potential.ab_initio.configuration_interaction.space != .cas) {
                std.log.err("TRUNCATED CI IS NOT SUPPORTED FOR NONADIABATIC DYNAMICS", .{});

                return error.InvalidInput;
            }
        }

        const is_cas = switch (opt.potential.ab_initio) {
            .casscf, .configuration_interaction => true,
            else => false,
        };

        if (is_cas) {
            const mult, const penalty = switch (opt.potential.ab_initio) {
                .casscf => |c| .{ c.multiplicity, c.spin_penalty },
                .configuration_interaction => |ci| .{ ci.space.cas.multiplicity, ci.space.cas.spin_penalty },
                else => unreachable,
            };

            if (mult > 1) {
                std.log.err("ONLY SINGLET MULTIPLICITY IS SUPPORTED FOR AB INITIO NAMD", .{});

                return error.InvalidInput;
            }

            if (penalty == null) {
                std.log.err("SPIN PENALTY IS MANDATORY FOR AB INITIO NAMD TO PREVENT TRIPLET CONTAMINATION", .{});

                return error.InvalidInput;
            }
        }

        switch (opt.nonadiabatic.?) {
            .ehrenfest => |eh| {
                const is_valid = switch (opt.potential.ab_initio) {
                    .casscf, .configuration_interaction => true,
                    else => false,
                };

                if (!is_valid) {
                    std.log.err("AB INITIO EHRENFEST IS ONLY SUPPORTED WITH CASSCF OR CONFIGURATION INTERACTION", .{});

                    return error.InvalidInput;
                }

                const nstate, const has_nacv = switch (opt.potential.ab_initio) {
                    .casscf => |c| .{ c.nstate, c.nacv != null },
                    .configuration_interaction => |ci| .{ ci.nstate, ci.nacv != null },
                    else => unreachable,
                };

                if (nstate < 2) {
                    std.log.err("AT LEAST 2 STATES ARE REQUIRED FOR EHRENFEST", .{});

                    return error.InvalidInput;
                }

                if (eh.coupling) |c| if (c != .nacv) {
                    std.log.err("ONLY NACV COUPLING IS SUPPORTED FOR AB INITIO EHRENFEST", .{});

                    return error.InvalidInput;
                };

                if (!has_nacv) {
                    std.log.err("NACV MUST BE CONFIGURED FOR EHRENFEST", .{});

                    return error.InvalidInput;
                }
            },

            .surface_hopping => |sh| switch (sh) {
                .fewest_switches => |fs| {
                    const is_valid = switch (opt.potential.ab_initio) {
                        .casscf, .configuration_interaction => true,
                        else => false,
                    };

                    if (!is_valid) {
                        std.log.err("AB INITIO FSSH IS ONLY SUPPORTED WITH CASSCF OR CONFIGURATION INTERACTION", .{});

                        return error.InvalidInput;
                    }

                    const nstate, const has_nacv = switch (opt.potential.ab_initio) {
                        .casscf => |c| .{ c.nstate, c.nacv != null },
                        .configuration_interaction => |ci| .{ ci.nstate, ci.nacv != null },
                        else => unreachable,
                    };

                    if (nstate < 2) {
                        std.log.err("AT LEAST 2 STATES ARE REQUIRED FOR FSSH", .{});

                        return error.InvalidInput;
                    }

                    if (fs.coupling) |c| if (c != .nacv) {
                        std.log.err("ONLY NACV COUPLING IS SUPPORTED FOR AB INITIO FSSH", .{});

                        return error.InvalidInput;
                    };

                    if (!has_nacv) {
                        std.log.err("NACV MUST BE CONFIGURED FOR FSSH", .{});

                        return error.InvalidInput;
                    }
                },

                .landau_zener => {
                    const is_valid = switch (opt.potential.ab_initio) {
                        .casscf, .configuration_interaction => true,
                        else => false,
                    };

                    if (!is_valid) {
                        std.log.err("LANDAU-ZENER ONLY SUPPORTS CASSCF OR CONFIGURATION INTERACTION", .{});

                        return error.InvalidInput;
                    }

                    const nstate = switch (opt.potential.ab_initio) {
                        .casscf => |c| c.nstate,
                        .configuration_interaction => |ci| ci.nstate,
                        else => unreachable,
                    };

                    if (nstate < 2) {
                        std.log.err("AT LEAST 2 STATES ARE REQUIRED FOR LANDAU-ZENER", .{});

                        return error.InvalidInput;
                    }
                },
            },
        }
    }

    if (opt.nonadiabatic) |na| switch (na) {
        .surface_hopping => |sh| switch (sh) {
            .fewest_switches => |fs| {
                if (fs.alpha != null and !opt.adiabatic) {
                    std.log.err("DECOHERENCE CORRECTION IS ONLY SUPPORTED FOR ADIABATIC SIMULATIONS", .{});

                    return error.InvalidInput;
                }

                if (fs.alpha) |a| if (std.math.isNan(a) or a <= 0) {
                    std.log.err("DECOHERENCE PARAMETER ALPHA MUST BE POSITIVE", .{});

                    return error.InvalidInput;
                };

                if (opt.adiabatic and fs.coupling == null) {
                    std.log.err("TIME-DERIVATIVE COUPLING MUST BE EXPLICITLY PROVIDED FOR ADIABATIC FSSH", .{});

                    return error.InvalidInput;
                }

                if (!opt.adiabatic and fs.coupling != null) {
                    std.log.err("TIME-DERIVATIVE COUPLING IS NOT SUPPORTED FOR DIABATIC SIMULATIONS", .{});

                    return error.InvalidInput;
                }

                if (fs.reflect_frustrated and (fs.coupling == null or fs.coupling.? != .nacv)) {
                    std.log.err("FRUSTRATED HOP REFLECTION REQUIRES NACV TIME-DERIVATIVE COUPLING", .{});

                    return error.InvalidInput;
                }

                if (fs.norm_tolerance) |tol| if (std.math.isNan(tol) or tol <= 0) {
                    std.log.err("NORM TOLERANCE MUST BE POSITIVE", .{});

                    return error.InvalidInput;
                };
            },

            inline else => {},
        },

        .ehrenfest => |eh| {
            if (opt.adiabatic and eh.coupling == null) {
                std.log.err("TIME-DERIVATIVE COUPLING MUST BE EXPLICITLY PROVIDED FOR ADIABATIC EHRENFEST", .{});

                return error.InvalidInput;
            }

            if (!opt.adiabatic and eh.coupling != null) {
                std.log.err("TIME-DERIVATIVE COUPLING IS NOT SUPPORTED FOR DIABATIC SIMULATIONS", .{});

                return error.InvalidInput;
            }

            if (eh.norm_tolerance) |tol| if (std.math.isNan(tol) or tol <= 0) {
                std.log.err("NORM TOLERANCE MUST BE POSITIVE", .{});

                return error.InvalidInput;
            };
        },
    };

    const has_state_coef = opt.write.state_coefficient != null or opt.stream.state_coefficient != null;

    if (has_state_coef) {
        const has_coefs = if (opt.nonadiabatic) |na| switch (na) {
            .surface_hopping => |sh| sh == .fewest_switches,
            .ehrenfest => true,
        } else false;

        if (!has_coefs) {
            std.log.err("STATE COEFFICIENT OUTPUT IS ONLY SUPPORTED WITH FEWEST SWITCHES OR EHRENFEST", .{});

            return error.InvalidInput;
        }
    }

    if (opt.potential == .ab_initio) switch (opt.potential.ab_initio) {
        inline .casscf, .configuration_interaction => |cas_or_ci| {
            if (opt.initial_conditions.state >= cas_or_ci.nstate) {
                std.log.err("REQUESTED ELECTRONIC STATE EXCEEDS NUMBER OF ELECTRONIC STATES", .{});

                return error.InvalidInput;
            }

            if (cas_or_ci.gradient) |grad| {
                const grad_states = switch (grad) {
                    inline else => |g| g.states,
                };

                if (grad_states != null) {
                    std.log.err("GRADIENT STATE MUST NOT BE SPECIFIED FOR AB INITIO POTENTIAL", .{});

                    return error.InvalidInput;
                }
            }

            if (cas_or_ci.nacv) |nacv| {
                const nacv_pairs = switch (nacv) {
                    inline else => |n| n.pairs,
                };

                if (nacv_pairs != null) {
                    std.log.err("NACV PAIRS MUST NOT BE SPECIFIED FOR AB INITIO POTENTIAL", .{});

                    return error.InvalidInput;
                }
            }

            if (cas_or_ci.nacv != null and cas_or_ci.nstate < 2) {
                std.log.err("AT LEAST 2 STATES ARE REQUIRED FOR NACV CALCULATION", .{});

                return error.InvalidInput;
            }
        },

        .hartree_fock => if (opt.initial_conditions.state > 0) {
            std.log.err("REQUESTED ELECTRONIC STATE IS NOT SUPPORTED FOR HARTREE-FOCK POTENTIAL", .{});

            return error.InvalidInput;
        },

        .moller_plesset => if (opt.initial_conditions.state > 0) {
            std.log.err("REQUESTED ELECTRONIC STATE IS NOT SUPPORTED FOR MOLLER-PLESSET POTENTIAL", .{});

            return error.InvalidInput;
        },
    };

    if (opt.thermostat != null and opt.nonadiabatic != null) {
        std.log.err("THERMOSTATS ARE NOT SUPPORTED FOR NONADIABATIC DYNAMICS", .{});

        return error.InvalidInput;
    }

    if (opt.thermostat) |topt| switch (topt) {
        .berendsen => |bopt| {
            if (bopt.temperature < 0) {
                std.log.err("TEMPERATURE MUST BE NON-NEGATIVE", .{});

                return error.InvalidInput;
            }

            if (bopt.tau <= 0) {
                std.log.err("COUPLING CONSTANT TAU MUST BE POSITIVE", .{});

                return error.InvalidInput;
            }
        },
        .langevin => |lopt| {
            if (lopt.temperature < 0) {
                std.log.err("TEMPERATURE MUST BE NON-NEGATIVE", .{});

                return error.InvalidInput;
            }

            if (lopt.gamma < 0) {
                std.log.err("FRICTION COEFFICIENT MUST BE NON-NEGATIVE", .{});

                return error.InvalidInput;
            }
        },
    };
}

/// Validates determinant excitation spaces, state counts, and derivative requests for configuration interaction.
fn checkConfigurationInteraction(opt: configuration_interaction.Options) !void {
    if (opt.hartree_fock.integral_direct) {
        std.log.err("INTEGRAL DIRECT CALCULATION IS NOT SUPPORTED FOR CONFIGURATION INTERACTION", .{});

        return error.InvalidInput;
    }

    try checkWriteCalculated(opt);

    switch (opt.space) {
        .cas => |c| {
            if (opt.hartree_fock.generalized) {
                std.log.err("GENERALIZED HARTREE-FOCK IS NOT SUPPORTED FOR CASCI", .{});

                return error.InvalidInput;
            }

            if (c.active_electrons == 0 or c.active_orbitals == 0) {
                std.log.err("CASCI ACTIVE ELECTRONS AND ACTIVE ORBITALS MUST BE GREATER THAN 0", .{});

                return error.InvalidInput;
            }

            if (c.multiplicity == 0) {
                std.log.err("CASCI MULTIPLICITY MUST BE GREATER THAN 0", .{});

                return error.InvalidInput;
            }

            if (c.spin_penalty) |pen| if (pen < 0 or std.math.isNan(pen)) {
                std.log.err("CASCI SPIN PENALTY MUST BE NON-NEGATIVE", .{});

                return error.InvalidInput;
            };
        },
        .excitations => |ex| {
            if (ex.len == 0) {
                std.log.err("CI EXCITATIONS LIST MUST NOT BE EMPTY", .{});

                return error.InvalidInput;
            }

            for (ex) |k| if (k == 0) {
                std.log.err("CI EXCITATION LEVEL MUST BE GREATER THAN 0", .{});

                return error.InvalidInput;
            };
        },
    }

    if (opt.nthreads == 0) {
        std.log.err("THREAD COUNT MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.nstate == 0) {
        std.log.err("NUMBER OF CI STATES MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    try checkMultiState("CI", opt);
}

/// Validates self-consistent field convergence parameters, spin multiplicity, and integral direct settings.
fn checkHartreeFock(opt: hartree_fock.Options) !void {
    if (opt.integral_direct) {
        if (opt.response != null) {
            std.log.err("CPHF RESPONSE PROPERTIES ARE NOT SUPPORTED FOR INTEGRAL DIRECT HARTREE-FOCK", .{});

            return error.InvalidInput;
        }

        if (opt.write.integrals.coulomb != null or opt.write.integrals.coulomb_d1 != null) {
            std.log.err("COULOMB INTEGRALS WRITE REQUESTED BUT INTEGRAL DIRECT HARTREE-FOCK IS ENABLED", .{});

            return error.InvalidInput;
        }
    }

    if (opt.nthreads == 0) {
        std.log.err("THREAD COUNT MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    try checkWriteCalculated(opt);

    if (opt.gradient == null or opt.gradient.? != .analytic) {
        if (opt.write.integrals.coulomb_d1 != null) {
            std.log.err("COULOMB DERIVATIVE WRITE REQUESTED BUT COULOMB DERIVATIVE IS NOT CALCULATED", .{});

            return error.InvalidInput;
        }

        if (opt.write.integrals.hamiltonian_d1 != null) {
            std.log.err("HAMILTONIAN DERIVATIVE WRITE REQUESTED BUT HAMILTONIAN DERIVATIVE IS NOT CALCULATED", .{});

            return error.InvalidInput;
        }

        if (opt.write.integrals.kinetic_d1 != null) {
            std.log.err("KINETIC DERIVATIVE WRITE REQUESTED BUT KINETIC DERIVATIVE IS NOT CALCULATED", .{});

            return error.InvalidInput;
        }

        if (opt.write.integrals.nuclear_d1 != null) {
            std.log.err("NUCLEAR DERIVATIVE WRITE REQUESTED BUT NUCLEAR DERIVATIVE IS NOT CALCULATED", .{});

            return error.InvalidInput;
        }

        if (opt.write.integrals.overlap_d1 != null) {
            std.log.err("OVERLAP DERIVATIVE WRITE REQUESTED BUT OVERLAP DERIVATIVE IS NOT CALCULATED", .{});

            return error.InvalidInput;
        }
    }

    if (opt.multiplicity == 0) {
        std.log.err("MULTIPLICITY MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.multiplicity > 1 and !opt.generalized) {
        std.log.err("OPEN SHELL IS NOT SUPPORTED UNLESS GENERALIZED HARTREE-FOCK IS ENABLED", .{});

        return error.OnlyClosedShellSupported;
    }

    if (opt.system.len == 0) {
        std.log.err("MOLECULAR SYSTEM XYZ PATH IS EMPTY", .{});

        return error.InvalidInput;
    }

    if (opt.basis.len == 0) {
        std.log.err("BASIS SET G94 PATH IS EMPTY", .{});

        return error.InvalidInput;
    }

    if (opt.iterations == 0) {
        std.log.err("SCF ITERATIONS MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.threshold <= 0) {
        std.log.err("SCF CONVERGENCE THRESHOLD MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.dft) |d| {
        if (opt.gradient) |g| if (g == .analytic) {
            std.log.err("ANALYTIC GRADIENTS ARE NOT SUPPORTED FOR DFT CALCULATIONS", .{});

            return error.InvalidInput;
        };

        if (opt.response != null) {
            std.log.err("CPKS RESPONSE PROPERTIES ARE NOT SUPPORTED FOR DFT CALCULATIONS", .{});

            return error.InvalidInput;
        }

        if (d.grid.radial == 0) {
            std.log.err("DFT RADIAL GRID POINTS MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        }

        const valid_angular_points = switch (d.grid.angular) {
            50, 74, 110, 302, 590, 974 => true,

            else => false,
        };

        if (!valid_angular_points) {
            std.log.err("DFT ANGULAR GRID SIZE IS NOT SUPPORTED. USE 50, 74, 110, 302, 590, OR 974", .{});

            return error.InvalidInput;
        }

        if (d.exchange_correlation != null and (d.exchange != null or d.correlation != null)) {
            std.log.err("CANNOT SPECIFY BOTH EXCHANGE_CORRELATION AND INDIVIDUAL FUNCTIONALS", .{});

            return error.InvalidInput;
        }

        if (d.exchange_correlation == null and d.exchange == null and d.correlation == null) {
            std.log.err("MUST SPECIFY EITHER EXCHANGE_CORRELATION OR EXCHANGE AND/OR CORRELATION FUNCTIONALS", .{});

            return error.InvalidInput;
        }
    }

    if (opt.response) |r| {
        if (r.iterations == 0) {
            std.log.err("CPHF RESPONSE ITERATIONS MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        }

        if (r.threshold <= 0) {
            std.log.err("CPHF RESPONSE THRESHOLD MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        }
    }
}

/// Validates atomic basis paths, molecular coordinates, and integral computation versus print requests.
fn checkMolecularIntegrals(opt: molecular_integrals.Options) !void {
    if (opt.system.len == 0) {
        std.log.err("MOLECULAR SYSTEM XYZ PATH IS EMPTY", .{});

        return error.InvalidInput;
    }

    if (opt.basis.len == 0) {
        std.log.err("BASIS SET G94 PATH IS EMPTY", .{});

        return error.InvalidInput;
    }

    if (opt.nthreads == 0) {
        std.log.err("THREAD COUNT MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    inline for (@typeInfo(@TypeOf(opt.calculate)).@"struct".field_names) |name| {
        const upper = comptime switch (name[0]) {
            'c' => "COULOMB",
            'h' => "HAMILTONIAN",
            'k' => "KINETIC",
            'n' => "NUCLEAR",
            'o' => "OVERLAP",
            else => unreachable,
        };
        const kind = comptime if (std.mem.endsWith(u8, name, "_d1")) "DERIVATIVE" else "INTEGRALS";

        const target_name = upper ++ " " ++ kind;

        const msg_print = target_name ++ " PRINT REQUESTED BUT " ++ target_name ++ " ARE NOT CALCULATED";
        const msg_write = target_name ++ " WRITE REQUESTED BUT " ++ target_name ++ " ARE NOT CALCULATED";

        if (@field(opt.log, name) and !@field(opt.calculate, name)) {
            std.log.err(msg_print, .{});

            return error.InvalidInput;
        }

        if (@field(opt.write, name) != null and !@field(opt.calculate, name)) {
            std.log.err(msg_write, .{});

            return error.InvalidInput;
        }
    }
}

/// Validates perturbation theory order, reference state compatibility, and correlation energy options.
fn checkMollerPlesset(opt: moller_plesset.Options) !void {
    if (opt.hartree_fock.integral_direct) {
        std.log.err("INTEGRAL DIRECT CALCULATION IS NOT SUPPORTED FOR MOLLER-PLESSET PERTURBATION THEORY", .{});

        return error.InvalidInput;
    }

    try checkWriteCalculated(opt);

    if (opt.order < 2) {
        std.log.err("MØLLER-PLESSET PERTURBATION ORDER MUST BE AT LEAST 2 (MP2)", .{});

        return error.InvalidInput;
    }

    if (opt.nthreads == 0) {
        std.log.err("THREAD COUNT MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }
}

/// Validates multi-state target index bounds, derivative state selections, and nonadiabatic coupling pairs.
fn checkMultiState(comptime method_name: []const u8, opt: anytype) !void {
    if (opt.gradient) |gradopt| {
        const states = switch (gradopt) {
            inline else => |g| g.states orelse &[_]u32{0},
        };

        if (states.len == 0) {
            std.log.err(method_name ++ " GRADIENT STATES LIST MUST NOT BE EMPTY", .{});

            return error.InvalidInput;
        }

        for (states) |s| if (s >= opt.nstate) {
            std.log.err("GRADIENT STATE MUST BE LESS THAN NUMBER OF STATES", .{});

            return error.InvalidInput;
        };
    }

    if (opt.hessian) |hessopt| {
        const states = switch (hessopt) {
            inline else => |h| h.states orelse &[_]u32{0},
        };

        if (states.len == 0) {
            std.log.err(method_name ++ " HESSIAN STATES LIST MUST NOT BE EMPTY", .{});

            return error.InvalidInput;
        }

        for (states) |s| if (s >= opt.nstate) {
            std.log.err("HESSIAN STATE MUST BE LESS THAN NUMBER OF STATES", .{});

            return error.InvalidInput;
        };
    }

    if (opt.nacv) |nacvopt| {
        if (opt.nstate < 2) {
            std.log.err("AT LEAST 2 STATES ARE REQUIRED FOR NACV CALCULATION", .{});

            return error.InvalidInput;
        }

        const pairs = switch (nacvopt) {
            inline else => |n| n.pairs orelse &[_][2]u32{.{ 0, 1 }},
        };

        if (pairs.len == 0) {
            std.log.err(method_name ++ " NACV PAIRS LIST MUST NOT BE EMPTY", .{});

            return error.InvalidInput;
        }

        for (pairs) |pair| {
            if (pair[0] == pair[1]) {
                std.log.err("NACV STATES MUST BE DISTINCT", .{});

                return error.InvalidInput;
            }

            if (pair[0] >= opt.nstate or pair[1] >= opt.nstate) {
                std.log.err("NACV STATES MUST BE LESS THAN NUMBER OF STATES", .{});

                return error.InvalidInput;
            }
        }
    }

    if (opt.optimize) |optopt| {
        const gradopt = switch (optopt) {
            .bfgs => |b| b.gradient,
            .steepest_descent => |s| s.gradient,
        };

        const states = switch (gradopt) {
            inline else => |g| g.states orelse &[_]u32{0},
        };

        if (states.len != 1) {
            std.log.err("EXACTLY ONE STATE MUST BE SPECIFIED FOR GEOMETRY OPTIMIZATION", .{});

            return error.InvalidInput;
        }

        if (states[0] >= opt.nstate) {
            std.log.err("OPTIMIZATION STATE MUST BE LESS THAN NUMBER OF STATES", .{});

            return error.InvalidInput;
        }
    }
}

/// Validates coordinate grid boundaries, mass arrays, boundary absorption, and wavepacket initial conditions.
fn checkQuantumDynamics(opt: quantum_dynamics.Options) !void {
    if (opt.partial_waves) |pw| {
        if (opt.grid.bounds.len != 1) {
            std.log.err("PARTIAL WAVES REQUIRE EXACTLY 1 GRID DIMENSION", .{});

            return error.InvalidInput;
        }

        if (pw.step == 0) {
            std.log.err("PARTIAL WAVES STEP MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        }

        if (pw.log_interval == 0) {
            std.log.err("PARTIAL WAVES LOG INTERVAL MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        }

        if (opt.flux_analysis == null) {
            std.log.err("PARTIAL WAVES REQUIRE FLUX ANALYSIS TO BE CONFIGURED", .{});

            return error.InvalidInput;
        }

        if (opt.flux_analysis) |fa| if (fa.flux_bounds.len > 0 and fa.flux_bounds[0][1] <= 0) {
            std.log.err("PARTIAL WAVES REQUIRE POSITIVE FLUX BOUND", .{});

            return error.InvalidInput;
        };
    }

    if (opt.potential == .ab_initio) {
        std.log.err("AB INITIO POTENTIAL IS NOT SUPPORTED FOR QUANTUM DYNAMICS", .{});

        return error.InvalidInput;
    }

    if (opt.time_step <= 0) {
        std.log.err("TIME STEP MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    for (opt.mass) |m| if (m <= 0) {
        std.log.err("MASS MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    };

    if (opt.mass.len != opt.initial_conditions.position.len) {
        std.log.err("MASS VECTOR MUST HAVE THE SAME LENGTH AS POSITION VECTOR", .{});

        return error.InvalidInput;
    }

    if (opt.log_interval == 0) {
        std.log.err("LOG INTERVAL MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.grid.npoint == 0) {
        std.log.err("GRID POINT COUNT MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.grid.bounds.len == 0) {
        std.log.err("GRID BOUNDS MUST NOT BE EMPTY", .{});

        return error.InvalidInput;
    }

    for (opt.grid.bounds) |b| if (b[0] >= b[1]) {
        std.log.err("GRID BOUNDS MIN MUST BE LESS THAN MAX", .{});

        return error.InvalidInput;
    };

    if (opt.initial_conditions.position.len == 0) {
        std.log.err("INITIAL POSITION VECTOR MUST NOT BE EMPTY", .{});

        return error.InvalidInput;
    }

    if (opt.initial_conditions.momentum.len != opt.initial_conditions.position.len) {
        std.log.err("INITIAL MOMENTUM AND POSITION VECTORS MUST HAVE THE SAME LENGTH", .{});

        return error.InvalidInput;
    }

    if (opt.initial_conditions.gamma.len != opt.initial_conditions.position.len) {
        std.log.err("INITIAL GAMMA VECTOR MUST HAVE THE SAME LENGTH AS POSITION VECTOR", .{});

        return error.InvalidInput;
    }

    if (opt.initial_conditions.position.len != opt.grid.bounds.len) {
        std.log.err("INITIAL CONDITIONS DIMENSION DOES NOT MATCH GRID DIMENSION", .{});

        return error.InvalidInput;
    }

    if (opt.absorbing_potential) |cap| {
        if (cap.bounds.len != opt.grid.bounds.len) {
            std.log.err("ABSORBING POTENTIAL BOUNDS DIMENSION DOES NOT MATCH GRID DIMENSION", .{});

            return error.InvalidInput;
        }

        for (cap.bounds, opt.grid.bounds) |cb, gb| if (cb[0] < gb[0] or cb[1] > gb[1] or cb[0] >= cb[1]) {
            std.log.err("ABSORBING POTENTIAL BOUNDS MUST LIE WITHIN GRID BOUNDS AND CB[0] < CB[1]", .{});

            return error.InvalidInput;
        };

        if (cap.exponent <= 0) {
            std.log.err("ABSORBING POTENTIAL EXPONENT MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        }

        if (cap.stop_norm) |sn| if (sn <= 0 or sn >= 1) {
            std.log.err("ABSORBING POTENTIAL STOP_NORM MUST BE BETWEEN 0 AND 1", .{});

            return error.InvalidInput;
        };
    }

    if (opt.imaginary) |imag| if (imag.nstate == 0) {
        std.log.err("ITP TARGET STATE COUNT MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    };

    if (opt.spectrum) |spec| if (spec.threshold <= 0 or spec.threshold >= 1) {
        std.log.err("SPECTRUM THRESHOLD MUST BE BETWEEN 0 AND 1", .{});

        return error.InvalidInput;
    };

    if (opt.flux_analysis) |flux| {
        for (opt.initial_conditions.momentum[1..]) |m| if (m != 0) {
            std.log.err("FLUX ANALYSIS REQUIRES ONLY THE FIRST MOMENTUM COORDINATE TO BE NON-ZERO", .{});

            return error.InvalidInput;
        };

        if (flux.flux_bounds.len != opt.grid.bounds.len) {
            std.log.err("FLUX ANALYSIS BOUNDS DIMENSION DOES NOT MATCH GRID DIMENSION", .{});

            return error.InvalidInput;
        }

        for (flux.flux_bounds, opt.grid.bounds) |fb, gb| if (fb[0] < gb[0] or fb[1] > gb[1] or fb[0] > fb[1]) {
            std.log.err("FLUX BOUNDS MUST LIE WITHIN GRID BOUNDS AND FB[0] <= FB[1]", .{});

            return error.InvalidInput;
        };

        if (flux.e_min < 0) {
            std.log.err("FLUX ANALYSIS E_MIN MUST BE NON-NEGATIVE", .{});

            return error.InvalidInput;
        }

        if (flux.e_max <= flux.e_min) {
            std.log.err("FLUX ANALYSIS E_MAX MUST BE GREATER THAN E_MIN", .{});

            return error.InvalidInput;
        }

        if (flux.asymptotic_offset >= flux.e_max) {
            std.log.err("FLUX ANALYSIS ASYMPTOTIC OFFSET MUST BE LESS THAN E_MAX", .{});

            return error.InvalidInput;
        }

        if (flux.npoint < 2) {
            std.log.err("FLUX ANALYSIS NPOINT MUST BE AT LEAST 2", .{});

            return error.InvalidInput;
        }

        if (flux.logarithmic and flux.e_min <= 0) {
            std.log.err("LOGARITHMIC FLUX ANALYSIS REQUIRES POSITIVE E_MIN", .{});

            return error.InvalidInput;
        }
    }

    if (opt.j_quantum_number != 0 and opt.grid.bounds.len != 1) {
        std.log.err("NONZERO J QUANTUM NUMBER REQUIRES EXACTLY 1 GRID DIMENSION", .{});

        return error.InvalidInput;
    }

    if (opt.grid.cylindrical) {
        if (opt.j_quantum_number != 0 or opt.partial_waves != null) {
            std.log.err("CYLINDRICAL SIMULATION DOES NOT SUPPORT NONZERO J QUANTUM NUMBER", .{});

            return error.InvalidInput;
        }

        if (opt.grid.bounds.len < 2) {
            std.log.err("CYLINDRICAL SIMULATION REQUIRES AT LEAST 2 GRID DIMENSIONS", .{});

            return error.InvalidInput;
        }

        const radial_dim = opt.grid.bounds.len - 1;

        if (opt.grid.bounds[radial_dim][0] != -opt.grid.bounds[radial_dim][1]) {
            std.log.err("CYLINDRICAL SIMULATION REQUIRES SYMMETRIC RADIAL GRID BOUNDS", .{});

            return error.InvalidInput;
        }

        if (opt.initial_conditions.position[radial_dim] != 0) {
            std.log.err("CYLINDRICAL SIMULATION REQUIRES RADIAL POSITION TO BE ZERO", .{});

            return error.InvalidInput;
        }

        if (opt.initial_conditions.momentum[radial_dim] != 0) {
            std.log.err("CYLINDRICAL SIMULATION REQUIRES RADIAL MOMENTUM TO BE ZERO", .{});

            return error.InvalidInput;
        }

        if (opt.absorbing_potential) |cap| {
            if (cap.bounds[radial_dim][0] != -cap.bounds[radial_dim][1]) {
                std.log.err("CYLINDRICAL SIMULATION REQUIRES SYMMETRIC RADIAL ABSORBING POTENTIAL BOUNDS", .{});

                return error.InvalidInput;
            }
        }

        if (opt.flux_analysis) |flux| {
            if (flux.flux_bounds[radial_dim][0] != -flux.flux_bounds[radial_dim][1]) {
                std.log.err("CYLINDRICAL SIMULATION REQUIRES SYMMETRIC RADIAL FLUX BOUNDS", .{});

                return error.InvalidInput;
            }
        }
    }
}

/// Verifies that disk export requests have corresponding gradient, Hessian, or nonadiabatic coupling targets.
fn checkWriteCalculated(opt: anytype) !void {
    if (opt.write.gradient != null and opt.gradient == null) {
        std.log.err("GRADIENT WRITE REQUESTED BUT GRADIENT IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.write.hessian != null and opt.hessian == null) {
        std.log.err("HESSIAN WRITE REQUESTED BUT HESSIAN IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (comptime @hasField(@TypeOf(opt.write), "nacv") and @hasField(@TypeOf(opt), "nacv")) {
        if (opt.write.nacv != null and opt.nacv == null) {
            std.log.err("NACV WRITE REQUESTED BUT NACV IS NOT CALCULATED", .{});

            return error.InvalidInput;
        }
    }

    if (comptime @hasField(@TypeOf(opt), "frequency") and @hasField(@TypeOf(opt), "hessian")) {
        if (opt.frequency != null and opt.hessian == null) {
            std.log.err("FREQUENCY CALCULATION REQUESTED BUT HESSIAN IS NOT CALCULATED", .{});

            return error.InvalidInput;
        }
    }
}
