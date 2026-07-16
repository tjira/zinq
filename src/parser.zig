//! Maps parameters from argument space $\mathcal{P}$ to physical Hamiltonian simulation targets.

const std = @import("std");

const Allocator = std.mem.Allocator;

const main = @import("main.zig");
const hartree_fock = @import("hartree_fock.zig");

const printf = @import("read_write.zig").printf;

/// Option representations for subcommands executed in physical basis space.
pub const Subcommand = enum {
    hf,
};

/// Represents subcommand arguments mapped to physical actions.
pub const SubcommandAction = struct {
    name: Subcommand,
    args: []const []const u8,
};

/// Represents the action mapping command parameters to the simulation trajectory.
pub const Action = union(enum) {
    help: void,
    files: []const []const u8,
    subcommand: SubcommandAction,
};

/// Maps input configuration paths and execution flags to the Hamiltonian state space.
pub const Parser = struct {
    action: Action,

    /// Help message for the main parser, detailing usage, options, and subcommands.
    pub const help_main =
        \\
        \\USAGE: zinq [INPUTS/SUBCOMMANDS] [OPTIONS]
        \\
        \\INPUTS:
        \\  file          JSON INPUT FILE DESCRIBING SIMULATION PARAMETERS (DEFAULT: input.json)
        \\
        \\OPTIONS:
        \\  -h, --help    PRINT THIS HELP MESSAGE AND EXIT
        \\
        \\SUBCOMMANDS:
        \\  hf            RUN HARTREE-FOCK METHOD DIRECTLY ON MOLECULAR COORDINATES
        \\
    ;

    /// Help message for the Hartree-Fock subcommand, detailing usage, options, and arguments.
    pub const help_hf =
        \\
        \\USAGE: zinq hf [ARGUMENTS] [OPTIONS]
        \\
        \\OPTIONS:
        \\  -b, --basis    SPECIFY BASIS SET (DEFAULT: STO-3G)
        \\  -h, --help     PRINT THIS HELP MESSAGE AND EXIT
        \\
        \\ARGUMENTS:
        \\  file           XYZ FILE DESCRIBING MOLECULE
        \\
    ;

    /// Initializes the parser by mapping the raw argument space $\mathcal{P}$ to a parsed action.
    pub fn init(args: []const []const u8, allocator: Allocator) !@This() {
        const action = try parse(args, allocator);

        return .{ .action = action };
    }

    /// Directs the parsed argument action into the corresponding state space execution path.
    pub fn dispatch(self: @This(), io: std.Io, gpa: Allocator, arena: Allocator) !void {
        switch (self.action) {
            .help => try runHelp(io, help_main),
            .subcommand => |command| try runSubcommand(io, gpa, arena, command),
            .files => |files| try runFiles(io, gpa, arena, files),
        }
    }

    /// Iteratively simulates physical systems described by the parsed json files.
    pub fn runFiles(io: std.Io, gpa: Allocator, arena: Allocator, files: []const []const u8) !void {
        for (files) |e| try main.run(f64, io, e, gpa, arena);
    }

    /// Outputs the configuration option spectrum to guide simulation setup.
    pub fn runHelp(io: std.Io, message: []const u8) !void {
        try printf(io, "{s}", .{message});
    }

    /// Enforces specific subcommand constraints on physical parameter evaluation.
    pub fn runSubcommand(io: std.Io, gpa: Allocator, arena: Allocator, sub: SubcommandAction) !void {
        switch (sub.name) {
            .hf => try runHartreeFock(io, gpa, arena, sub),
        }
    }

    /// Projects CLI subcommand parameters into Hartree-Fock electronic states.
    pub fn runHartreeFock(io: std.Io, gpa: Allocator, arena: Allocator, sub: SubcommandAction) !void {
        var molecule: ?[]const u8 = null;
        var basis: []const u8 = "sto-3g";

        var i: usize = 0;

        while (i < sub.args.len) : (i += 1) {
            const arg = sub.args[i];

            if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
                try runHelp(io, help_hf);

                return;
            }

            if (std.mem.eql(u8, arg, "-b") or std.mem.eql(u8, arg, "--basis")) {
                if (i + 1 >= sub.args.len) {
                    try printf(io, "MISSING VALUE FOR BASIS OPTION\n", .{});

                    return error.MissingBasisValue;
                }

                basis, i = .{ sub.args[i + 1], i + 1 };

                continue;
            }

            if (std.mem.startsWith(u8, arg, "-")) {
                try printf(io, "UNKNOWN '{s}' OPTION\n", .{arg});

                return error.UnknownOption;
            }

            if (molecule != null) {
                try printf(io, "MULTIPLE MOLECULE FILES SPECIFIED\n", .{});

                return error.MultipleMoleculeFiles;
            }

            molecule = arg;
        }

        if (molecule == null) {
            try printf(io, "MOLECULE FILE IS REQUIRED FOR 'hf' SUBCOMMAND\n", .{});

            return error.MissingMoleculeFile;
        }

        const basis_resolved = try std.fmt.allocPrint(arena, "builtin:{s}", .{basis});

        const opt = hartree_fock.Options{ .system = molecule.?, .basis = basis_resolved };

        var result = try hartree_fock.run(f64, io, opt, true, gpa);
        defer result.deinit(gpa);
    }

    /// Projects the raw command line token sequence into distinct execution pathways.
    fn parse(args: []const []const u8, allocator: Allocator) !Action {
        if (args.len <= 1) {
            const default_files = try allocator.alloc([]const u8, 1);

            default_files[0] = "input.json";

            return .{ .files = default_files };
        }

        if (std.mem.eql(u8, args[1], "hf")) {
            return .{ .subcommand = .{ .name = .hf, .args = args[2..] } };
        }

        for (args[1..]) |arg| if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            return .help;
        };

        var files: std.ArrayList([]const u8) = .empty;
        errdefer files.deinit(allocator);

        for (args[1..]) |arg| {
            try files.append(allocator, arg);
        }

        return .{ .files = try files.toOwnedSlice(allocator) };
    }
};
