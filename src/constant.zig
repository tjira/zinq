//! PHYSICAL CONSTANTS ACCORDING TO CODATA 2022
const std = @import("std");

pub const c = 299792458.0; // SPEED OF LIGHT [m/s]
pub const a0 = 5.29177210544e-11; // BOHR RADIUS [m]
pub const Eh = 4.3597447222060e-18; // HARTREE ENERGY [J]
pub const mu = 1.66053906892e-27; // ATOMIC MASS UNIT [kg]

pub const A2BOHR = 1e-10 / a0; // ANGSTROM TO BOHR
pub const AU2CM = @sqrt(Eh / (mu * a0 * a0)) / (2e2 * c * std.math.pi); // ATOMIC UNIT OF FREQUENCY TO CM^-1

pub const AN2SM = std.StaticStringMap(i32).initComptime(.{
    .{ "H", 1 },
    .{ "He", 2 },
    .{ "Li", 3 },
    .{ "Be", 4 },
    .{ "B", 5 },
    .{ "C", 6 },
    .{ "N", 7 },
    .{ "O", 8 },
    .{ "F", 9 },
    .{ "Ne", 10 },
    .{ "Na", 11 },
    .{ "Mg", 12 },
    .{ "Al", 13 },
    .{ "Si", 14 },
    .{ "P", 15 },
    .{ "S", 16 },
    .{ "Cl", 17 },
    .{ "Ar", 18 },
});

pub const AN2M = std.StaticStringMap(f64).initComptime(.{
    .{ "H", 1.007825032 },
    .{ "He", 4.002602 },
    .{ "Li", 7.016004 },
    .{ "Be", 9.012182 },
    .{ "B", 11.009305 },
    .{ "C", 12.000000 },
    .{ "N", 14.003074 },
    .{ "O", 15.994915 },
    .{ "F", 18.998403 },
    .{ "Ne", 19.99244 },
    .{ "Na", 22.989770 },
    .{ "Mg", 23.98504 },
    .{ "Al", 26.981538 },
    .{ "Si", 27.976927 },
    .{ "P", 30.973762 },
    .{ "S", 31.972071 },
    .{ "Cl", 34.968853 },
    .{ "Ar", 39.962383 },
});
