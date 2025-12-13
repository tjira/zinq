//! Error context for capturing the first error in a parallel execution environment.

const std = @import("std");

pub const ErrorContext = struct {
    mutex: std.Thread.Mutex = .{},
    err: ?anyerror = null,

    pub fn capture(self: *@This(), err: anyerror) void {
        self.mutex.lock(); defer self.mutex.unlock();

        if (self.err == null) self.err = err;
    }
};
