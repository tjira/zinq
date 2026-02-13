//! Expression evaluator target

const std = @import("std");

const device_write = @import("device_write.zig");
const error_handling = @import("error_handling.zig");
const shunting_yard = @import("shunting_yard.zig");

const print = device_write.print;
const printJson = device_write.printJson;
const shuntingYard = shunting_yard.shuntingYard;
const throw = error_handling.throw;

/// Options for the expression evaluator target.
pub fn Options(comptime _: type) type {
    return struct {
        expression: []const u8,
        variables: []const []const u8,
        values: []const f64
    };
}

/// Output structure for the expression evaluator target.
pub fn Output(comptime T: type) type {
    return struct {
        value: T,

        /// Free the output structure.
        pub fn deinit(_: @This(), _: std.mem.Allocator) void {}
    };
}

/// Run the expression evaluator target with the provided options and allocator, returning the output structure.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    if (opt.variables.len != opt.values.len) return throw(Output(T), "VARIABLES AND VALUES LENGTH MISMATCH", .{});

    if (enable_printing) try print("\nEXPRESSION: '{s}'", .{opt.expression});

    if (enable_printing and opt.variables.len > 0) {

        try print(", VARIABLES: ", .{});

        for (0..opt.variables.len) |i| try print("{s}={d}{s}", .{opt.variables[i], opt.values[i], if (i == opt.variables.len - 1) "" else ", "});
    }

    var rpn = try shuntingYard(T, opt.expression, opt.variables, allocator); defer rpn.deinit(allocator);

    var map = std.StringHashMap(T).init(allocator); defer map.deinit();

    for (opt.variables, opt.values) |variable, value| try map.put(variable, value);

    const rpnstr = try rpn.toString(allocator); defer allocator.free(rpnstr); const value = try rpn.evaluate(map);

    if (enable_printing) try print("\n\nRPN: {s}\n\nVALUE: {d:.8}\n", .{rpnstr, value});

    return Output(T){.value = value};
}
