const std = @import("std");

pub fn Vector(comptime T: type) type {
    return struct {
        data: []T,
        shape: [1]usize,

        pub fn init(size: usize, gpa: std.mem.Allocator) !@This() {
            return @This(){
                .data = try gpa.alloc(T, size),
                .shape = .{size},
            };
        }

        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }
    };
}
