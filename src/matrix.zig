pub fn Matrix(comptime T: type) type {
    return struct {
        data: []T,
        shape: [2]usize,
    };
}
