pub fn Vector(comptime T: type) type {
    return struct {
        data: []T,
        shape: [1]usize,
    };
}
