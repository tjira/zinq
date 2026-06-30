//! PHYSICAL CONSTANTS ACCORDING TO CODATA 2022
const std = @import("std");

pub const a0 = 5.29177210544e-11; // BOHR RADIUS [m]

pub const A2BOHR = 1e-10 / a0; // ANGSTROM TO BOHR

pub const AN2SM = std.StaticStringMap(i32).initComptime(.{
    // zig fmt: off
    .{ "H",   1 },
    .{ "He",  2 },
    .{ "Li",  3 },
    .{ "Be",  4 },
    .{ "B",   5 },
    .{ "C",   6 },
    .{ "N",   7 },
    .{ "O",   8 },
    .{ "F",   9 },
    .{ "Ne", 10 },
    .{ "Na", 11 },
    .{ "Mg", 12 },
    .{ "Al", 13 },
    .{ "Si", 14 },
    .{ "P",  15 },
    .{ "S",  16 },
    .{ "Cl", 17 },
    .{ "Ar", 18 },
    // zig fmt: on
});
