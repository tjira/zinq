//! Basis library.

const std = @import("std");

const contracted_gaussian = @import("contracted_gaussian.zig");
const classical_particle = @import("classical_particle.zig");
const embedded_files = @import("embedded_files.zig");
const global_variables = @import("global_variables.zig");

const ContractedGaussian = contracted_gaussian.ContractedGaussian;
const ClassicalParticle = classical_particle.ClassicalParticle;

const AN2SM = global_variables.AN2SM;
const BASIS_FILES = embedded_files.BASIS_FILES;

/// Basis struct.
pub fn BasisSet(comptime T: type) type {
    return struct {
        contracted_gaussians: []const ContractedGaussian(T), array: std.ArrayList(T), name: []const u8,

        /// Get the basis set for the given system and name.
        pub fn init(system: ClassicalParticle(T), name: []const u8, allocator: std.mem.Allocator) !BasisSet(T) {
            var basis = std.ArrayList(ContractedGaussian(T)){}; errdefer basis.deinit(allocator);

            var array = std.ArrayList(T){}; errdefer array.deinit(allocator);

            const basis_file_contents = BASIS_FILES.get(name) orelse {

                std.log.err("BASIS SET '{s}' NOT FOUND", .{name});

                return error.InvalidInput;
            };

            const basis_json = try std.json.parseFromSlice(std.json.Value, allocator, basis_file_contents, .{}); defer basis_json.deinit();

            for (0..system.atoms.?.len) |i| {

                var atomic_number_string: [3]u8 = undefined;

                const atomic_number_length = (try std.fmt.bufPrint(&atomic_number_string, "{d}", .{system.atoms.?[i]})).len;

                if (basis_json.value.object.get("elements").?.object.get(atomic_number_string[0..atomic_number_length]) == null) {

                    std.log.err("BASIS SET '{s}' NOT DEFINED FOR ATOMIC NUMBER {d}", .{name, system.atoms.?[i]});

                    return error.InvalidInput;
                }

                const shells = basis_json.value.object.get("elements").?.object.get(atomic_number_string[0..atomic_number_length]).?.object.get("electron_shells").?.array.items;

                for (shells) |shell| {

                    const am        = shell.object.get("angular_momentum").?.array.items[0].integer;
                    const exponents = shell.object.get("exponents"       ).?.array.items;
                    const coefs     = shell.object.get("coefficients"    ).?.array.items;

                    var a     = try allocator.alloc(T, 3            ); defer allocator.free(a    );
                    var c     = try allocator.alloc(T, exponents.len); defer allocator.free(c    );
                    var alpha = try allocator.alloc(T, exponents.len); defer allocator.free(alpha);

                    for (coefs) |coef| {

                        for (0..alpha.len) |k| c[k]     = try std.fmt.parseFloat(T, coef.array.items[k].string);
                        for (0..alpha.len) |k| alpha[k] = try std.fmt.parseFloat(T, exponents       [k].string);

                        try array.append(allocator, @as(T, @floatFromInt(alpha.len)));
                        try array.append(allocator, @as(T, @floatFromInt(am       )));
                        
                        try array.append(allocator, system.position.data[3 * i + 0]);
                        try array.append(allocator, system.position.data[3 * i + 1]);
                        try array.append(allocator, system.position.data[3 * i + 2]);

                        for (alpha) |val| try array.append(allocator, val);
                        for (c    ) |val| try array.append(allocator, val);

                        for (0..@as(usize, @intCast(am + 1))) |j| {

                            const lx = @as(usize, @intCast(am)) - j;

                            for (0..@as(usize, @intCast(am + 1)) - lx) |k| {

                                const ly = @as(usize, @intCast(am)) - lx - k;

                                a[0] = @as(T, @floatFromInt(lx));
                                a[1] = @as(T, @floatFromInt(ly));

                                a[2] = @as(T, @floatFromInt(am)) - @as(T, @floatFromInt(lx + ly));

                                const center = .{system.position.data[3 * i], system.position.data[3 * i + 1], system.position.data[3 * i + 2]};

                                const an = try std.fmt.parseInt(usize, atomic_number_string[0..atomic_number_length], 10);

                                try basis.append(allocator, try ContractedGaussian(T).init(an, center, .{a[0], a[1], a[2]}, c, alpha, allocator));
                            }
                        }
                    }
                }
            }

            return .{.contracted_gaussians = try basis.toOwnedSlice(allocator), .name = name, .array = array};
        }

        /// Deinitialize the basis set.
        pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
            for (self.contracted_gaussians) |cg| cg.deinit(allocator);

            self.array.deinit(allocator);

            allocator.free(self.contracted_gaussians);
        }

        /// Get the number of basis functions.
        pub fn nbf(self: BasisSet(T)) usize {
            return self.contracted_gaussians.len;
        }
    };
}
