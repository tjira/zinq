//! This module provides functions for calculating integrals using the libint library.

const std = @import("std");

const config = @import("config");

const libint = if (config.use_libint) @cImport(@cInclude("libint.h")) else struct {};

const basis_set = @import("basis_set.zig");
const classical_particle = @import("classical_particle.zig");
const real_matrix = @import("real_matrix.zig");
const real_tensor_four = @import("real_tensor_four.zig");

const BasisSet = basis_set.BasisSet;
const ClassicalParticle = classical_particle.ClassicalParticle;
const RealMatrix = real_matrix.RealMatrix;
const RealTensor4 = real_tensor_four.RealTensor4;

/// Calculate the Coulomb integrals.
pub fn coulomb(comptime T: type, I: *RealTensor4(T), system: ClassicalParticle(T), basis: BasisSet(T)) !void {
    if (comptime !config.use_libint) @compileError("LIBINT IS NOT ENABLED");

    libint.coulomb(&I.data[0], system.atoms.?.len, &system.atoms.?[0], &system.position.data[0], basis.array.items.len, basis.array.items.ptr);
}

/// Calculate the kinetic integrals.
pub fn kinetic(comptime T: type, I: *RealMatrix(T), system: ClassicalParticle(T), basis: BasisSet(T)) !void {
    if (comptime !config.use_libint) @compileError("LIBINT IS NOT ENABLED");

    libint.kinetic(&I.data[0], system.atoms.?.len, &system.atoms.?[0], &system.position.data[0], basis.array.items.len, basis.array.items.ptr);
}

/// Calculate the nuclear integrals.
pub fn nuclear(comptime T: type, I: *RealMatrix(T), system: ClassicalParticle(T), basis: BasisSet(T)) !void {
    if (comptime !config.use_libint) @compileError("LIBINT IS NOT ENABLED");

    libint.nuclear(&I.data[0], system.atoms.?.len, &system.atoms.?[0], &system.position.data[0], basis.array.items.len, basis.array.items.ptr);
}

/// Calculate the overlap integrals.
pub fn overlap(comptime T: type, I: *RealMatrix(T), system: ClassicalParticle(T), basis: BasisSet(T)) !void {
    if (comptime !config.use_libint) @compileError("LIBINT IS NOT ENABLED");

    libint.overlap(&I.data[0], system.atoms.?.len, &system.atoms.?[0], &system.position.data[0], basis.array.items.len, basis.array.items.ptr);
}
