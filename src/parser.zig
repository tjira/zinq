//! Maps parameters from argument space $\mathcal{P}$ to physical Hamiltonian simulation targets.

const std = @import("std");

const Allocator = std.mem.Allocator;

const main = @import("main.zig");
const hartree_fock = @import("hartree_fock.zig");
const moller_plesset = @import("moller_plesset.zig");

const printf = @import("read_write.zig").printf;

/// Represents the action mapping command parameters to the simulation trajectory.
pub const Action = union(enum) {
    help: void,
    files: []const []const u8,
    subcommand: SubcommandAction,
};

/// Represents the phase space mapping of discrete physical options and coordinate arguments.
pub const ParsedArgs = struct {
    positional: []const []const u8,
    options: std.StringHashMap([]const u8),
};

/// Option representations for subcommands executed in physical basis space.
pub const Subcommand = enum {
    hf,
    mp,
};

/// Represents subcommand arguments mapped to physical actions.
pub const SubcommandAction = struct {
    name: Subcommand,
    args: []const []const u8,
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
        \\  mp            RUN MOLLER-PLESSET PERTURBATION THEORY ON MOLECULAR COORDINATES
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

    /// Help message for the Moller-Plesset subcommand, detailing usage, options, and arguments.
    pub const help_mp =
        \\
        \\USAGE: zinq mp [ARGUMENTS] [OPTIONS]
        \\
        \\OPTIONS:
        \\  -b, --basis    SPECIFY BASIS SET (DEFAULT: STO-3G)
        \\  -o, --order    SPECIFY PERTURBATION ORDER (DEFAULT: 2)
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
            .mp => try runMollerPlesset(io, gpa, arena, sub),
        }
    }

    /// Projects CLI subcommand parameters into Hartree-Fock electronic states.
    pub fn runHartreeFock(io: std.Io, gpa: Allocator, arena: Allocator, sub: SubcommandAction) !void {
        const parsed = try parseArgs(sub.args, arena);

        if (parsed.options.contains("-h") or parsed.options.contains("--help")) {
            try runHelp(io, help_hf);

            return;
        }

        var opt_iter = parsed.options.keyIterator();

        while (opt_iter.next()) |key| if (!std.mem.eql(u8, key.*, "-b") and !std.mem.eql(u8, key.*, "--basis")) {
            try printf(io, "UNKNOWN '{s}' OPTION\n", .{key.*});

            return error.UnknownOption;
        };

        if (parsed.positional.len > 1) {
            try printf(io, "MULTIPLE MOLECULE FILES SPECIFIED\n", .{});

            return error.MultipleMoleculeFiles;
        }

        if (parsed.positional.len == 0) {
            try printf(io, "MOLECULE FILE IS REQUIRED FOR 'hf' SUBCOMMAND\n", .{});

            return error.MissingMoleculeFile;
        }

        const basis_opt = parsed.options.get("-b") orelse parsed.options.get("--basis");

        if (basis_opt) |b| if (b.len == 0) {
            try printf(io, "MISSING VALUE FOR BASIS OPTION\n", .{});

            return error.MissingBasisValue;
        };

        const basis_resolved = try std.fmt.allocPrint(arena, "builtin:{s}", .{basis_opt orelse "sto-3g"});

        const opt = hartree_fock.Options{
            .system = parsed.positional[0],
            .basis = basis_resolved,
        };

        var result = try hartree_fock.run(f64, io, opt, true, gpa);
        defer result.deinit(gpa);
    }

    /// Projects Hartree-Fock reference states into perturbed Møller-Plesset correlation spaces.
    pub fn runMollerPlesset(io: std.Io, gpa: Allocator, arena: Allocator, sub: SubcommandAction) !void {
        const parsed = try parseArgs(sub.args, arena);

        if (parsed.options.contains("-h") or parsed.options.contains("--help")) {
            try runHelp(io, help_mp);

            return;
        }

        var opt_iter = parsed.options.keyIterator();

        while (opt_iter.next()) |key| {
            const is_basis = std.mem.eql(u8, key.*, "-b") or std.mem.eql(u8, key.*, "--basis");
            const is_order = std.mem.eql(u8, key.*, "-o") or std.mem.eql(u8, key.*, "--order");

            if (!is_basis and !is_order) {
                try printf(io, "UNKNOWN '{s}' OPTION\n", .{key.*});

                return error.UnknownOption;
            }
        }

        if (parsed.positional.len > 1) {
            try printf(io, "MULTIPLE MOLECULE FILES SPECIFIED\n", .{});

            return error.MultipleMoleculeFiles;
        }

        if (parsed.positional.len == 0) {
            try printf(io, "MOLECULE FILE IS REQUIRED FOR 'mp' SUBCOMMAND\n", .{});

            return error.MissingMoleculeFile;
        }

        const basis_opt = parsed.options.get("-b") orelse parsed.options.get("--basis");

        if (basis_opt) |b| if (b.len == 0) {
            try printf(io, "MISSING VALUE FOR BASIS OPTION\n", .{});

            return error.MissingBasisValue;
        };

        const order_opt = parsed.options.get("-o") orelse parsed.options.get("--order");

        if (order_opt) |o| if (o.len == 0) {
            try printf(io, "MISSING VALUE FOR ORDER OPTION\n", .{});

            return error.MissingOrderValue;
        };

        const order = std.fmt.parseInt(u32, order_opt orelse "2", 10) catch |err| {
            try printf(io, "INVALID VALUE FOR ORDER OPTION\n", .{});

            return err;
        };

        const basis_resolved = try std.fmt.allocPrint(arena, "builtin:{s}", .{basis_opt orelse "sto-3g"});

        const opt = moller_plesset.Options{
            .hartree_fock = .{
                .system = parsed.positional[0],
                .basis = basis_resolved,
            },
            .order = order,
        };

        var result = try moller_plesset.run(f64, io, opt, true, gpa);
        defer result.deinit(gpa);
    }

    /// Projects coordinate-space trajectories and physical configuration options from raw token streams.
    fn parseArgs(args: []const []const u8, allocator: Allocator) !ParsedArgs {
        var positional: std.ArrayList([]const u8), var options = .{ .empty, std.StringHashMap([]const u8).init(allocator) };

        var i: usize = 0;

        while (i < args.len) : (i += 1) {
            const arg = args[i];

            if (std.mem.startsWith(u8, arg, "-")) {
                if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
                    try options.put(arg, "");

                    continue;
                }

                if (i + 1 >= args.len) try options.put(arg, "");

                if (i + 1 < args.len) {
                    try options.put(arg, args[i + 1]);

                    i += 1;
                }

                continue;
            }

            try positional.append(allocator, arg);
        }

        return .{ .positional = try positional.toOwnedSlice(allocator), .options = options };
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

        if (std.mem.eql(u8, args[1], "mp")) {
            return .{ .subcommand = .{ .name = .mp, .args = args[2..] } };
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
