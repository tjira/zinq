//! Functions related to parallel processing and concurrency.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");

const throw = error_handling.throw;

const PARALLEL_ERROR = &global_variables.PARALLEL_ERROR;

/// Checks if a parallel error has occurred and throws it if so.
pub fn checkParallelError() !void {
    if (PARALLEL_ERROR.* != null) try throw(void, "ERROR IN PARALLEL EXECUTION, CHECK STDERR FOR DETAILS", .{});
}
