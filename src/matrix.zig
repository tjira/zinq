const std = @import("std");

pub fn Matrix(comptime T: type) type {
    return struct {
        data: []T,
        shape: [2]usize,

        pub fn init(rows: usize, cols: usize, gpa: std.mem.Allocator) !@This() {
            return @This(){
                .data = try gpa.alloc(T, rows * cols),
                .shape = .{ rows, cols },
            };
        }

        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }
    };
}
