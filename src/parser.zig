//! Maps parameters from argument space $\mathcal{P}$ to physical Hamiltonian simulation targets.

const std = @import("std");

const Allocator = std.mem.Allocator;
const ArrayList = std.ArrayList;
const StringHashMap = std.StringHashMap;

const main = @import("main.zig");
const hartree_fock = @import("hartree_fock.zig");
const moller_plesset = @import("moller_plesset.zig");
const tensor = @import("tensor.zig");

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
        \\  mm            PERFORM MATRIX MULTIPLICATION ON INPUT MATRIX FILES
        \\  mp            RUN MOLLER-PLESSET PERTURBATION THEORY ON MOLECULAR COORDINATES
        \\  randn         GENERATE PSEUDORANDOM GAUSSIAN MATRICES
        \\
    ;

    /// Help message for the Hartree-Fock subcommand, detailing usage, options, and arguments.
    pub const help_hf =
        \\
        \\USAGE: zinq hf [ARGUMENTS] [OPTIONS]
        \\
        \\OPTIONS:
        \\  -b, --basis           SPECIFY BASIS SET (DEFAULT: STO-3G)
        \\  -s, --multiplicity    SPECIFY MULTIPLICITY (DEFAULT: 1)
        \\  -c, --charge          SPECIFY CHARGE (DEFAULT: 0)
        \\  --generalized         ENABLE GENERALIZED HARTREE-FOCK
        \\  -h, --help            PRINT THIS HELP MESSAGE AND EXIT
        \\
        \\ARGUMENTS:
        \\  file                  XYZ FILE DESCRIBING MOLECULE
        \\
    ;

    /// Help message for the matrix multiplication subcommand, detailing usage, options, and arguments.
    pub const help_mm =
        \\
        \\USAGE: zinq mm [ARGUMENTS] [OPTIONS]
        \\
        \\OPTIONS:
        \\  -o, --output          OUTPUT FILE PATH TO SAVE RESULT MATRIX
        \\  -a, --alpha           SCALAR MULTIPLIER ALPHA (DEFAULT: 1)
        \\  --trans-a             TRANSPOSE FIRST MATRIX A
        \\  --trans-b             TRANSPOSE SECOND MATRIX B
        \\  -p, --print           PRINT RESULT MATRIX TO TERMINAL
        \\  -h, --help            PRINT THIS HELP MESSAGE AND EXIT
        \\
        \\ARGUMENTS:
        \\  file_a                FIRST INPUT MATRIX FILE
        \\  file_b                SECOND INPUT MATRIX FILE
        \\
    ;

    /// Help message for the Moller-Plesset subcommand, detailing usage, options, and arguments.
    pub const help_mp =
        \\
        \\USAGE: zinq mp [ARGUMENTS] [OPTIONS]
        \\
        \\OPTIONS:
        \\  -b, --basis           SPECIFY BASIS SET (DEFAULT: STO-3G)
        \\  -o, --order           SPECIFY PERTURBATION ORDER (DEFAULT: 2)
        \\  -s, --multiplicity    SPECIFY MULTIPLICITY (DEFAULT: 1)
        \\  -c, --charge          SPECIFY CHARGE (DEFAULT: 0)
        \\  --generalized         ENABLE GENERALIZED PERTURBATION THEORY
        \\  -h, --help            PRINT THIS HELP MESSAGE AND EXIT
        \\
        \\ARGUMENTS:
        \\  file                  XYZ FILE DESCRIBING MOLECULE
        \\
    ;

    /// Help message for standard normal random matrix generator subcommand, detailing options and arguments.
    pub const help_randn =
        \\
        \\USAGE: zinq randn [ARGUMENTS] [OPTIONS]
        \\
        \\OPTIONS:
        \\  -o, --output          OUTPUT FILE PATH TO SAVE GENERATED MATRIX
        \\  -s, --seed            SPECIFY PSEUDORANDOM SEED (DEFAULT: 0)
        \\  -p, --print           PRINT GENERATED MATRIX TO TERMINAL
        \\  -h, --help            PRINT THIS HELP MESSAGE AND EXIT
        \\
        \\ARGUMENTS:
        \\  rows                  NUMBER OF MATRIX ROWS
        \\  cols                  NUMBER OF MATRIX COLUMNS
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

    /// Projects CLI subcommand parameters into Hartree-Fock electronic states.
    pub fn runHartreeFock(io: std.Io, gpa: Allocator, arena: Allocator, sub: SubcommandAction) !void {
        const allowed_options = &.{
            "-b",
            "--basis",
            "-c",
            "--charge",
            "-s",
            "--multiplicity",
        };

        const allowed_flags = &.{
            "--generalized",
        };

        const parsed = try parseArgs(io, sub.args, arena, allowed_options, allowed_flags);

        if (parsed.options.contains("-h") or parsed.options.contains("--help")) {
            try runHelp(io, help_hf);

            return;
        }

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

        const multiplicity_opt = parsed.options.get("-s") orelse parsed.options.get("--multiplicity");

        if (multiplicity_opt) |s| if (s.len == 0) {
            try printf(io, "MISSING VALUE FOR MULTIPLICITY OPTION\n", .{});

            return error.MissingMultiplicityValue;
        };

        const charge_opt = parsed.options.get("-c") orelse parsed.options.get("--charge");

        if (charge_opt) |c| if (c.len == 0) {
            try printf(io, "MISSING VALUE FOR CHARGE OPTION\n", .{});

            return error.MissingChargeValue;
        };

        const multiplicity = if (multiplicity_opt) |s| std.fmt.parseInt(u32, s, 10) catch |err| {
            try printf(io, "INVALID VALUE FOR MULTIPLICITY OPTION\n", .{});

            return err;
        } else 1;

        const charge = if (charge_opt) |c| std.fmt.parseInt(i32, c, 10) catch |err| {
            try printf(io, "INVALID VALUE FOR CHARGE OPTION\n", .{});

            return err;
        } else 0;

        const basis_resolved = try std.fmt.allocPrint(arena, "builtin:{s}", .{basis_opt orelse "sto-3g"});

        const opt = hartree_fock.Options{
            .system = parsed.positional[0],
            .basis = basis_resolved,
            .multiplicity = multiplicity,
            .charge = charge,
            .generalized = parsed.options.contains("--generalized"),
        };

        var result = try hartree_fock.run(f64, io, opt, true, gpa);
        defer result.deinit(gpa);
    }

    /// Outputs the configuration option spectrum to guide simulation setup.
    pub fn runHelp(io: std.Io, message: []const u8) !void {
        try printf(io, "{s}", .{message});
    }

    /// Executes matrix multiplication $\mathbf{C} = \alpha \mathbf{A} \mathbf{B}$ from parsed command line options.
    pub fn runMatmul(io: std.Io, gpa: Allocator, arena: Allocator, sub: SubcommandAction) !void {
        const allowed_options = &.{
            "-a",
            "--alpha",
            "-o",
            "--output",
        };

        const allowed_flags = &.{
            "-p",
            "--print",
            "--trans-a",
            "--trans-b",
        };

        const parsed = try parseArgs(io, sub.args, arena, allowed_options, allowed_flags);

        if (parsed.options.contains("-h") or parsed.options.contains("--help")) {
            try runHelp(io, help_mm);

            return;
        }

        if (parsed.positional.len > 2) {
            try printf(io, "MULTIPLE MATRIX FILES SPECIFIED\n", .{});

            return error.MultipleMatrixFiles;
        }

        if (parsed.positional.len < 2) {
            try printf(io, "TWO MATRIX FILES ARE REQUIRED FOR 'mm' SUBCOMMAND\n", .{});

            return error.MissingMatrixFile;
        }

        const alpha_opt = parsed.options.get("-a") orelse parsed.options.get("--alpha");

        if (alpha_opt) |a| if (a.len == 0) {
            try printf(io, "MISSING VALUE FOR ALPHA OPTION\n", .{});

            return error.MissingAlphaValue;
        };

        const output_opt = parsed.options.get("-o") orelse parsed.options.get("--output");

        if (output_opt) |o| if (o.len == 0) {
            try printf(io, "MISSING VALUE FOR OUTPUT OPTION\n", .{});

            return error.MissingOutputValue;
        };

        const alpha = if (alpha_opt) |a| std.fmt.parseFloat(f64, a) catch |err| {
            try printf(io, "INVALID VALUE FOR ALPHA OPTION\n", .{});

            return err;
        } else 1;

        const opt = tensor.MatmulOptions{
            .a = parsed.positional[0],
            .b = parsed.positional[1],
            .alpha = alpha,
            .log = .{ .product = parsed.options.contains("-p") or parsed.options.contains("--print") },
            .trans_a = parsed.options.contains("--trans-a"),
            .trans_b = parsed.options.contains("--trans-b"),
            .write = .{ .product = output_opt },
        };

        var result = try tensor.runMatmul(f64, io, opt, true, gpa);
        defer result.deinit(gpa);
    }

    /// Projects Hartree-Fock reference states into perturbed Møller-Plesset correlation spaces.
    pub fn runMollerPlesset(io: std.Io, gpa: Allocator, arena: Allocator, sub: SubcommandAction) !void {
        const allowed_options = &.{
            "-b",
            "--basis",
            "-c",
            "--charge",
            "-o",
            "--order",
            "-s",
            "--multiplicity",
        };

        const allowed_flags = &.{
            "--generalized",
        };

        const parsed = try parseArgs(io, sub.args, arena, allowed_options, allowed_flags);

        if (parsed.options.contains("-h") or parsed.options.contains("--help")) {
            try runHelp(io, help_mp);

            return;
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

        const multiplicity_opt = parsed.options.get("-s") orelse parsed.options.get("--multiplicity");

        if (multiplicity_opt) |s| if (s.len == 0) {
            try printf(io, "MISSING VALUE FOR MULTIPLICITY OPTION\n", .{});

            return error.MissingMultiplicityValue;
        };

        const charge_opt = parsed.options.get("-c") orelse parsed.options.get("--charge");

        if (charge_opt) |c| if (c.len == 0) {
            try printf(io, "MISSING VALUE FOR CHARGE OPTION\n", .{});

            return error.MissingChargeValue;
        };

        const multiplicity = if (multiplicity_opt) |s| std.fmt.parseInt(u32, s, 10) catch |err| {
            try printf(io, "INVALID VALUE FOR MULTIPLICITY OPTION\n", .{});

            return err;
        } else 1;

        const charge = if (charge_opt) |c| std.fmt.parseInt(i32, c, 10) catch |err| {
            try printf(io, "INVALID VALUE FOR CHARGE OPTION\n", .{});

            return err;
        } else 0;

        const basis_resolved = try std.fmt.allocPrint(arena, "builtin:{s}", .{basis_opt orelse "sto-3g"});

        const opt = moller_plesset.Options{
            .hartree_fock = .{
                .system = parsed.positional[0],
                .basis = basis_resolved,
                .multiplicity = multiplicity,
                .charge = charge,
                .generalized = parsed.options.contains("--generalized"),
            },
            .order = order,
        };

        var result = try moller_plesset.run(f64, io, opt, true, gpa);
        defer result.deinit(gpa);
    }

    /// Generates pseudorandom tensor elements sampled from normal distribution $\mathcal{N}(0, 1)$.
    pub fn runRandn(io: std.Io, gpa: Allocator, arena: Allocator, sub: SubcommandAction) !void {
        const allowed_options = &.{
            "-o",
            "--output",
            "-s",
            "--seed",
        };

        const allowed_flags = &.{
            "-p",
            "--print",
        };

        const parsed = try parseArgs(io, sub.args, arena, allowed_options, allowed_flags);

        if (parsed.options.contains("-h") or parsed.options.contains("--help")) {
            try runHelp(io, help_randn);

            return;
        }

        if (parsed.positional.len > 2) {
            try printf(io, "TOO MANY ARGUMENTS SPECIFIED FOR 'randn'\n", .{});

            return error.TooManyArguments;
        }

        if (parsed.positional.len < 2) {
            try printf(io, "SHAPE (ROWS AND COLS) IS REQUIRED FOR 'randn' SUBCOMMAND\n", .{});

            return error.MissingShapeArgument;
        }

        const output_opt = parsed.options.get("-o") orelse parsed.options.get("--output");

        if (output_opt) |o| if (o.len == 0) {
            try printf(io, "MISSING VALUE FOR OUTPUT OPTION\n", .{});

            return error.MissingOutputValue;
        };

        const seed_opt = parsed.options.get("-s") orelse parsed.options.get("--seed");

        if (seed_opt) |s| if (s.len == 0) {
            try printf(io, "MISSING VALUE FOR SEED OPTION\n", .{});

            return error.MissingSeedValue;
        };

        const rows = std.fmt.parseInt(usize, parsed.positional[0], 10) catch |err| {
            try printf(io, "INVALID VALUE FOR ROWS ARGUMENT\n", .{});

            return err;
        };

        const cols = std.fmt.parseInt(usize, parsed.positional[1], 10) catch |err| {
            try printf(io, "INVALID VALUE FOR COLS ARGUMENT\n", .{});

            return err;
        };

        const seed = if (seed_opt) |s| std.fmt.parseInt(u64, s, 10) catch |err| {
            try printf(io, "INVALID VALUE FOR SEED OPTION\n", .{});

            return err;
        } else 0;

        const opt = tensor.RandomOptions{
            .distribution = .{ .normal = .{} },
            .log = .{ .matrix = parsed.options.contains("-p") or parsed.options.contains("--print") },
            .seed = seed,
            .shape = .{ rows, cols },
            .write = .{ .matrix = output_opt },
        };

        var result = try tensor.runRandom(f64, io, opt, true, gpa);
        defer result.deinit(gpa);
    }

    /// Enforces specific subcommand constraints on physical parameter evaluation.
    pub fn runSubcommand(io: std.Io, gpa: Allocator, arena: Allocator, sub: SubcommandAction) !void {
        switch (sub.name) {
            .hf => try runHartreeFock(io, gpa, arena, sub),
            .mm => try runMatmul(io, gpa, arena, sub),
            .mp => try runMollerPlesset(io, gpa, arena, sub),
            .randn => try runRandn(io, gpa, arena, sub),
        }
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

        if (std.mem.eql(u8, args[1], "mm")) {
            return .{ .subcommand = .{ .name = .mm, .args = args[2..] } };
        }

        if (std.mem.eql(u8, args[1], "mp")) {
            return .{ .subcommand = .{ .name = .mp, .args = args[2..] } };
        }

        if (std.mem.eql(u8, args[1], "randn")) {
            return .{ .subcommand = .{ .name = .randn, .args = args[2..] } };
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

    /// Maps trajectories and options from raw token streams to allowed options in space $\mathcal{P}$.
    fn parseArgs(io: std.Io, args: []const []const u8, allocator: Allocator, comptime options_allowed: []const []const u8, comptime flags_allowed: []const []const u8) !ParsedArgs {
        var positional, var options = .{ ArrayList([]const u8).empty, StringHashMap([]const u8).init(allocator) };

        var i: usize = 0;

        while (i < args.len) : (i += 1) {
            const arg = args[i];

            if (std.mem.startsWith(u8, arg, "-")) {
                if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
                    try options.put(arg, "");

                    continue;
                }

                var is_flag = false;

                inline for (flags_allowed) |opt| if (std.mem.eql(u8, arg, opt)) {
                    is_flag = true;

                    break;
                };

                if (is_flag) {
                    try options.put(arg, "");

                    continue;
                }

                var is_option = false;

                inline for (options_allowed) |opt| if (std.mem.eql(u8, arg, opt)) {
                    is_option = true;

                    break;
                };

                if (!is_option) {
                    try printf(io, "UNKNOWN '{s}' OPTION\n", .{arg});

                    return error.UnknownOption;
                }

                if (i + 1 >= args.len or std.mem.startsWith(u8, args[i + 1], "-")) {
                    try options.put(arg, "");
                }

                if (i + 1 < args.len and !std.mem.startsWith(u8, args[i + 1], "-")) {
                    try options.put(arg, args[i + 1]);

                    i += 1;
                }

                continue;
            }

            try positional.append(allocator, arg);
        }

        return .{ .positional = try positional.toOwnedSlice(allocator), .options = options };
    }
};

/// Option representations for subcommands executed in physical basis space.
pub const Subcommand = enum {
    hf,
    mm,
    mp,
    randn,
};

/// Represents subcommand arguments mapped to physical actions.
pub const SubcommandAction = struct {
    name: Subcommand,
    args: []const []const u8,
};
