//! File that contains parallel algorithms for Zig.

const std = @import("std");

const error_handling = @import("error_handling.zig");

const throw = error_handling.throw;

/// Function for parallel execution of a for loop.
pub fn parallelFor(comptime func: anytype, args: anytype, start: usize, end: usize, nthread: usize, allocator: std.mem.Allocator) !void {
    if (nthread == 0) return throw(void, "NUMBER OF THREADS MUST BE GREATER THAN ZERO", .{});
    if (end <= start) return throw(void, "END MUST BE GREATER THAN START IN FOR LOOP", .{});

    var pool: std.Thread.Pool = undefined; try pool.init(.{.n_jobs = nthread, .allocator = allocator}); defer pool.deinit();

    for (start..end) |i| try pool.spawn(func, .{i, args});
}
