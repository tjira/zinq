//! Basis library.

const std = @import("std");

const contracted_gaussian = @import("contracted_gaussian.zig");
const classical_particle = @import("classical_particle.zig");
const embedded_files = @import("embedded_files.zig");
const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");

const ContractedGaussian = contracted_gaussian.ContractedGaussian;
const ClassicalParticle = classical_particle.ClassicalParticle;

const throw = error_handling.throw;

const AN2SM = global_variables.AN2SM;
const BASIS_FILES = embedded_files.BASIS_FILES;

/// Basis struct.
pub fn BasisSet(comptime T: type) type {
    return struct {
        contracted_gaussians: []const ContractedGaussian(T),

        allocator: std.mem.Allocator,

        /// Get the basis set for the given system and name.
        pub fn init(system: ClassicalParticle(T), name: []const u8, allocator: std.mem.Allocator) !BasisSet(T) {
            var basis = std.ArrayList(ContractedGaussian(T)){};

            const basis_file_contents = BASIS_FILES.get(name) orelse return throw(BasisSet(T), "BASIS SET \"{s}\" NOT FOUND", .{name});

            const basis_json = try std.json.parseFromSlice(std.json.Value, allocator, basis_file_contents, .{}); defer basis_json.deinit();

            for (0..system.atoms.?.len) |i| {

                var atomic_number_string: [3]u8 = undefined;

                const atomic_number_length = (try std.fmt.bufPrint(&atomic_number_string, "{d}", .{system.atoms.?[i]})).len;

                if (basis_json.value.object.get("elements").?.object.get(atomic_number_string[0..atomic_number_length]) == null) {
                    return throw(BasisSet(T), "ATOMIC NUMBER {d} ({s}) NOT FOUND IN BASIS SET FILE", .{system.atoms.?[i], try AN2SM(system.atoms.?[i])});
                }

                const shells = basis_json.value.object.get("elements").?.object.get(atomic_number_string[0..atomic_number_length]).?.object.get("electron_shells").?.array.items;

                for (shells) |shell| {

                    const am = shell.object.get("angular_momentum").?.array.items[0].integer;
                    const exponents = shell.object.get("exponents").?.array.items;
                    const coefs = shell.object.get("coefficients").?.array.items;

                    var a = try allocator.alloc(T, 3); defer allocator.free(a);
                    var c = try allocator.alloc(T, exponents.len); defer allocator.free(c);
                    var alpha = try allocator.alloc(T, exponents.len); defer allocator.free(alpha);

                    for (coefs) |coef| {

                        for (0..alpha.len) |k| c[k] = try std.fmt.parseFloat(T, coef.array.items[k].string);

                        for (0..@as(usize, @intCast(am + 1))) |lx| {
                            for (0..@as(usize, @intCast(am + 1)) - lx) |ly| {

                                a[0] = @as(T, @floatFromInt(lx));
                                a[1] = @as(T, @floatFromInt(ly));
                                a[2] = @as(T, @floatFromInt(am)) - @as(T, @floatFromInt(lx + ly));

                                for (0..alpha.len) |k| alpha[k] = try std.fmt.parseFloat(T, exponents[k].string);

                                const center = .{system.position.data[3 * i], system.position.data[3 * i + 1], system.position.data[3 * i + 2]};

                                try basis.append(allocator, try ContractedGaussian(T).init(center, .{a[0], a[1], a[2]}, c, alpha, allocator));
                            }
                        }
                    }
                }
            }

            return .{.contracted_gaussians = try basis.toOwnedSlice(allocator), .allocator = allocator};
        }

        /// Deinitialize the basis set.
        pub fn deinit(self: BasisSet(T)) void {
            for (self.contracted_gaussians) |cg| cg.deinit();

            self.allocator.free(self.contracted_gaussians);
        }

        /// Get the number of basis functions.
        pub fn nbf(self: BasisSet(T)) usize {
            return self.contracted_gaussians.len;
        }
    };
}
