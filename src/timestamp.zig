//! File that provides a function to get the current timestamp.

const std = @import("std");

const global_variables = @import("global_variables.zig");

const EPOCH_YEARS_OFFSET = global_variables.EPOCH_YEARS_OFFSET;

/// Get the current timestamp in the ISO 8601 format.
pub fn timestamp(allocator: std.mem.Allocator) ![]const u8 {
    const mts: u64 = @intCast(std.time.milliTimestamp());

    const year = getYear(mts);
    const month = getMonth(mts);
    const day = getDay(mts);
    const hour = getHour(mts);
    const minute = getMinute(mts);
    const second = getSecond(mts);
    const milli = getMilli(mts);

    return try std.fmt.allocPrint(allocator, "{d:04}-{d:02}-{d:02}T{d:02}:{d:02}:{d:02}.{d:03}Z", .{year, month, day, hour, minute, second, milli});
}

/// Get the day (1-31) from milliseconds timestamp since epoch.
pub fn getDay(mts: u64) u5 {
    var days = @divTrunc(mts, std.time.ms_per_day); var year: std.time.epoch.Year = EPOCH_YEARS_OFFSET;

    while (days >= std.time.epoch.getDaysInYear(year)) {
        days -= std.time.epoch.getDaysInYear(year); year += 1;
    }

    var month: std.time.epoch.Month = .jan;

    while (days >= std.time.epoch.getDaysInMonth(year, month)) {
        days -= std.time.epoch.getDaysInMonth(year, month); month = @enumFromInt(@intFromEnum(month) + 1);
    }

    return @intCast(days + 1);
}

/// Get the hour (0-23) from milliseconds timestamp since epoch.
pub fn getHour(mts: u64) u5 {
    return @intCast(@divTrunc(mts, std.time.ms_per_hour) % 24);
}

/// Get the millisecond (0-999) from milliseconds timestamp since epoch.
pub fn getMilli(mts: u64) u10 {
    return @intCast(mts % std.time.ms_per_s);
}

/// Get the minute (0-59) from milliseconds timestamp since epoch.
pub fn getMinute(mts: u64) u6 {
    return @intCast(@divTrunc(mts, std.time.ms_per_min) % 60);
}

/// Get the month (1-12) from milliseconds timestamp since epoch.
pub fn getMonth(mts: u64) u4 {
    var days = @divTrunc(mts, std.time.ms_per_day); var year: std.time.epoch.Year = EPOCH_YEARS_OFFSET;

    while (days >= std.time.epoch.getDaysInYear(year)) {
        days -= std.time.epoch.getDaysInYear(year); year += 1;
    }

    var month: std.time.epoch.Month = .jan;

    while (days >= std.time.epoch.getDaysInMonth(year, month)) {
        days -= std.time.epoch.getDaysInMonth(year, month); month = @enumFromInt(@intFromEnum(month) + 1);
    }

    return month.numeric();
}

/// Get the second (0-59) from milliseconds timestamp since epoch.
pub fn getSecond(mts: u64) u6 {
    return @intCast(@divTrunc(mts, std.time.ms_per_s) % 60);
}

/// Get the year from milliseconds timestamp since epoch.
pub fn getYear(mts: u64) std.time.epoch.Year {
    var days = @divTrunc(mts, std.time.ms_per_day); var year: std.time.epoch.Year = 1970;

    while (days >= std.time.epoch.getDaysInYear(year)) {
        days -= std.time.epoch.getDaysInYear(year); year += 1;
    }

    return year;
}
