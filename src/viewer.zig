//! Interactive graphical rendering interface for projected molecular geometry and orbital density manifolds.

const config = @import("config");
const std = @import("std");

const viewer = @import("cimport.zig").viewer;

/// Configuration options for the interactive coordinate visualization viewport.
pub const Options = struct {};

/// Initializes GLFW, binds the OpenGL context, and displays the empty viewport window until closed.
pub fn openWindow() bool {
    if (config.has_graphics) {
        return viewer.viewer_open_window();
    }

    return false;
}
