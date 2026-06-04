const std = @import("std");

pub fn div(a: anytype, b: anytype) @TypeOf(a, b) {
    return if (comptime @typeInfo(@TypeOf(a)) == .@"struct") a.div(b) else a / b;
}

pub fn mul(a: anytype, b: anytype) @TypeOf(a, b) {
    return if (comptime @typeInfo(@TypeOf(a)) == .@"struct") a.mul(b) else a * b;
}

pub fn is_struct(a: anytype) bool {
    return comptime @typeInfo(@TypeOf(a)) == .@"struct";
}
