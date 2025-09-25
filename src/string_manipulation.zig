//! Several functions for manipulationg strings.

pub fn uncr(string: []const u8) []const u8 {
    return if (string[string.len - 1] == 13) string[0..string.len - 1] else string;
}
