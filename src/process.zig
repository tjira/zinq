//! This module provides functionality to execute terminal commands and capture their output.

const std = @import("std");

/// Executes a provided string as a command in the terminal and returns its output.
pub fn executeCommand(io: std.Io, command: []const u8, allocator: std.mem.Allocator) ![]u8 {
    if (std.mem.containsAtLeast(u8, command, 1, "&&")) {
        var outputs: std.ArrayList(u8) = .empty;
        errdefer outputs.deinit(allocator);

        var iterator = std.mem.tokenizeAny(u8, command, "&&");

        while (iterator.next()) |cmd| {
            const output = try executeCommand(io, cmd, allocator);
            try outputs.appendSlice(allocator, output);
        }

        return try outputs.toOwnedSlice(allocator);
    }

    var args: std.ArrayList([]u8) = .empty;
    defer args.deinit(allocator);

    var iterator = std.mem.tokenizeAny(u8, command, " ");
    var string = false;

    while (iterator.next()) |token| {
        if (string) {
            const arg = try allocator.realloc(args.items[args.items.len - 1], args.items[args.items.len - 1].len + token.len + 1);

            arg[args.items[args.items.len - 1].len] = ' ';

            @memcpy(arg[args.items[args.items.len - 1].len + 1 ..], token);

            args.items[args.items.len - 1] = arg;
        } else try args.append(allocator, try allocator.dupe(u8, token));

        if (token[0] == '\'' or token[token.len - 1] == '\'') {
            var arg = args.items[args.items.len - 1];

            var start: usize = 0;
            var end: usize = arg.len;

            if (token[0] == '\'') start += 1;
            if (token[token.len - 1] == '\'') end -= 1;

            if (start > 0) std.mem.copyForwards(u8, arg[0 .. end - start], arg[start..end]);

            args.items[args.items.len - 1] = try allocator.realloc(arg, end - start);
        }

        if (token[0] == '\'') string = true;
        if (token[token.len - 1] == '\'') string = false;
    }

    const result = try std.process.run(allocator, io, .{.argv = args.items});

    for (args.items) |arg| allocator.free(arg);

    allocator.free(result.stderr);

    if (result.term.exited == 0) return result.stdout;

    return error.CommandExecutionFailed;
}
