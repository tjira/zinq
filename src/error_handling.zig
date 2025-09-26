//! Functions related to the error handling mechanism in Zig.

const std = @import("std");

const device_write = @import("device_write.zig");

const write = device_write.write;

/// Throws a runtime error with the given message.
pub fn throw(comptime T: type, comptime message: []const u8, args: anytype) !T {
    return throwSpecific(T, message, args, error.RuntimeError);
}

/// Throws a specific error with the given message.
pub fn throwSpecific(comptime T: type, comptime message: []const u8, args: anytype, err: anyerror) !T {

    try write(std.fs.File.stderr(), message, args);
    try write(std.fs.File.stderr(), " -- ", .{});

    return err;
}
