const std = @import("std");

const basis_names = [_][]const u8{
    "3-21g",
    "6-21g",
    "6-311g",
    "6-311gs",
    "6-311gss",
    "6-311pg",
    "6-311pgs",
    "6-311pgss",
    "6-311ppg",
    "6-311ppgs",
    "6-311ppgss",
    "6-31g",
    "6-31gs",
    "6-31gss",
    "6-31pg",
    "6-31pgs",
    "6-31pgss",
    "6-31ppg",
    "6-31ppgs",
    "6-31ppgss",
    "aug-cc-pv5z",
    "aug-cc-pvdz",
    "aug-cc-pvqz",
    "aug-cc-pvtz",
    "cc-pv5z",
    "cc-pvdz",
    "cc-pvqz",
    "cc-pvtz",
    "def2-qzvp",
    "def2-qzvpd",
    "def2-qzvpp",
    "def2-qzvppd",
    "def2-svp",
    "def2-svpd",
    "def2-tzvp",
    "def2-tzvpd",
    "def2-tzvpp",
    "def2-tzvppd",
    "sto-2g",
    "sto-3g",
    "sto-4g",
    "sto-5g",
    "sto-6g",
};

pub const bases = std.StaticStringMap([]const u8).initComptime(buildBases());

fn buildBases() [basis_names.len]struct { []const u8, []const u8 } {
    var result: [basis_names.len]struct { []const u8, []const u8 } = undefined;

    inline for (basis_names, 0..) |name, i| {
        result[i] = .{ name, @embedFile("basis/" ++ name ++ ".g94") };
    }

    return result;
}
