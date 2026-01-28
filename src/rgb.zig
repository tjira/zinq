//! RGB class implementation in Zig.

const std = @import("std");

/// RGB color representation.
pub const RGB = struct {
    r: u8, g: u8, b: u8
};
