//! Interactive graphical rendering interface for projected molecular geometry and orbital density manifolds.

const config = @import("config");
const std = @import("std");

const viewer = @import("cimport.zig").viewer;

/// Configuration options for the interactive coordinate visualization viewport.
pub const Options = struct {};

/// Initializes GLFW, binds the OpenGL context, and displays the molecular viewport window until closed.
pub fn openWindow(filename: []const u8, allocator: std.mem.Allocator) !bool {
    if (config.has_graphics) {
        const path_c = try allocator.dupeSentinel(u8, filename, 0);
        defer allocator.free(path_c);

        return viewer.viewer_open_window(path_c.ptr);
    }

    return false;
}
