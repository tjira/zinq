//! This module provides functionality to execute terminal commands and capture their output.

const std = @import("std");

const global_variables = @import("global_variables.zig");

const MAX_COMMAND_OUTPUT_BYTES = global_variables.MAX_COMMAND_OUTPUT_BYTES;

/// Executes a provided string as a command in the terminal and returns its output.
pub fn executeCommand(command: []const u8, allocator: std.mem.Allocator) ![]u8 {
    if (std.mem.containsAtLeast(u8, command, 1, "&&")) {

        var outputs = std.ArrayList(u8){}; errdefer outputs.deinit(allocator);

        var iterator = std.mem.tokenizeAny(u8, command, "&&");

        while (iterator.next()) |cmd| {
            const output = try executeCommand(cmd, allocator); try outputs.appendSlice(allocator, output);
        }

        return try outputs.toOwnedSlice(allocator);
    }

    var args = std.ArrayList([]u8){}; defer args.deinit(allocator);

    var iterator = std.mem.tokenizeAny(u8, command, " "); var string = false;

    while (iterator.next()) |token| {

        if (string) {

            const arg = try allocator.realloc(args.items[args.items.len - 1], args.items[args.items.len - 1].len + token.len + 1);
            
            arg[args.items[args.items.len - 1].len] = ' ';
            
            @memcpy(arg[args.items[args.items.len - 1].len + 1..], token);

            args.items[args.items.len - 1] = arg;

        } else try args.append(allocator, try allocator.dupe(u8, token));

        if (token[0] == '\'' or token[token.len - 1] == '\'') {

            var arg = args.items[args.items.len - 1];

            var start: usize = 0; var end: usize = arg.len;

            if (token[0] == '\'') start += 1; if (token[token.len - 1] == '\'') end -= 1;

            if (start > 0) std.mem.copyForwards(u8, arg[0..end - start], arg[start..end]);

            args.items[args.items.len - 1] = try allocator.realloc(arg, end - start);
        }

        if (token[0] == '\'') string = true; if (token[token.len - 1] == '\'') string = false;
    }

    var stdout = std.ArrayList(u8){}; var stderr = std.ArrayList(u8){}; errdefer stdout.deinit(allocator); defer stderr.deinit(allocator);

    var child = std.process.Child.init(args.items, allocator);

    child.stdout_behavior = .Pipe;
    child.stderr_behavior = .Pipe;

    try child.spawn(); try child.collectOutput(allocator, &stdout, &stderr, MAX_COMMAND_OUTPUT_BYTES);

    const term = try child.wait();

    for (args.items) |arg| allocator.free(arg);

    if (term.Exited == 0) return stdout.toOwnedSlice(allocator);

    return error.CommandExecutionFailed;
}
