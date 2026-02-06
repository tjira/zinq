//! This module provides functionality to execute terminal commands and capture their output.

const std = @import("std");

/// Executes a provided string as a command in the terminal and returns its output.
pub fn executeCommand(command: []const u8, allocator: std.mem.Allocator) ![]u8 {
    var args = std.ArrayList([]const u8){}; defer args.deinit(allocator);

    var iterator = std.mem.tokenizeAny(u8, command, " ");

    while (iterator.next()) |token| {
        try args.append(allocator, token);
    }

    var stdout = std.ArrayList(u8){}; var stderr = std.ArrayList(u8){}; defer stderr.deinit(allocator);

    var child = std.process.Child.init(args.items, allocator);

    child.stdout_behavior = .Pipe;
    child.stderr_behavior = .Pipe;

    try child.spawn(); try child.collectOutput(allocator, &stdout, &stderr, 1048576);

    const term = try child.wait();

    if (term.Exited == 0) return stdout.toOwnedSlice(allocator);

    return error.CommandExecutionFailed;
}
