const std = @import("std");

pub fn print(io: std.Io, str: []const u8) !void {
    try std.Io.File.stdout().writeStreamingAll(io, str);
}

pub fn printf(io: std.Io, comptime format: []const u8, args: anytype) !void {
    var buffer: [4096]u8 = undefined;
    var writer = std.Io.File.stdout().writer(io, &buffer);

    try writer.interface.print(format, args);
    try writer.interface.flush();
}
