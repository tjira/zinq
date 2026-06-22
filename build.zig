const std = @import("std");

pub fn build(b: *std.Build) !void {
    const target = b.standardTargetOptions(.{});

    const zinq_module = try setupZinq(b, b.standardOptimizeOption(.{}), target);

    setupTests(b, zinq_module);
}

fn setupZinq(b: *std.Build, opt: std.builtin.OptimizeMode, target: std.Build.ResolvedTarget) !*std.Build.Module {
    const zinq_module = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = opt,
        .strip = opt == .ReleaseFast or opt == .ReleaseSmall,
        .link_libc = true,
        .link_libcpp = true,
    });

    try linkDependencies(b, zinq_module);

    const exe_zinq = b.addExecutable(.{
        .name = "zinq",
        .root_module = zinq_module,
    });

    if (target.query.isNative()) {
        b.installArtifact(exe_zinq);
    }

    if (!target.query.isNative()) {
        const dest: std.Build.InstallDir = .{ .custom = try getTriple(b, target) };

        b.getInstallStep().dependOn(&b.addInstallArtifact(exe_zinq, .{ .dest_dir = .{ .override = dest } }).step);
    }

    const run_exe_zinq = b.addRunArtifact(exe_zinq);

    b.step("run", "Run the application").dependOn(&run_exe_zinq.step);

    return zinq_module;
}

fn setupTests(b: *std.Build, zinq_module: *std.Build.Module) void {
    const test_module = b.createModule(.{
        .root_source_file = b.path("test/main.zig"),
        .target = zinq_module.resolved_target,
        .strip = zinq_module.strip,
        .optimize = zinq_module.optimize,
        .link_libc = zinq_module.link_libc,
        .link_libcpp = true,
    });

    test_module.addImport("zinq", zinq_module);

    const exe_test = b.addTest(.{
        .root_module = test_module,
    });

    const run_exe_test = b.addRunArtifact(exe_test);

    b.step("test", "Run unit tests").dependOn(&run_exe_test.step);
}

fn linkDependencies(b: *std.Build, module: *std.Build.Module) !void {
    const dirs = [_][]const u8{ "lib", "include", "include/eigen3" };

    const triple, var archos: []const u8 = .{ try getTriple(b, module.resolved_target.?), undefined };

    if (std.mem.lastIndexOfScalar(u8, triple, '-')) |i| {
        archos = triple[0..i];
    }

    const ext = try std.fmt.allocPrint(b.allocator, "external-{s}", .{archos});

    std.Io.Dir.cwd().access(b.graph.io, ext, .{}) catch |err| {
        if (err == error.FileNotFound) {
            std.log.err("REQUIRED DEPENDENCIES DIRECTORY '{s}' DOES NOT EXIST", .{ext});
        }

        return err;
    };

    const dir0 = try std.fmt.allocPrint(b.allocator, "{s}/{s}", .{ ext, dirs[0] });
    const dir1 = try std.fmt.allocPrint(b.allocator, "{s}/{s}", .{ ext, dirs[1] });
    const dir2 = try std.fmt.allocPrint(b.allocator, "{s}/{s}", .{ ext, dirs[2] });

    module.addLibraryPath(.{ .cwd_relative = dir0 });
    module.addIncludePath(.{ .cwd_relative = dir1 });
    module.addIncludePath(.{ .cwd_relative = dir2 });

    module.addCSourceFile(.{ .file = b.path("src/libint.cpp") });

    const libint_translate = b.addTranslateC(.{
        .root_source_file = b.path("src/libint.h"),
        .target = module.resolved_target.?,
        .optimize = module.optimize.?,
    });

    libint_translate.addIncludePath(.{ .cwd_relative = dir1 });

    module.addImport("libint", libint_translate.createModule());

    const libs = [_][]const u8{
        "fftw3",
        "int2",
        "openblas",
        "xc",
    };

    for (libs) |lib| {
        module.linkSystemLibrary(lib, .{});
    }
}

fn getTriple(b: *std.Build, target: std.Build.ResolvedTarget) ![]const u8 {
    const triple = .{ @tagName(target.result.cpu.arch), @tagName(target.result.os.tag), @tagName(target.result.abi) };

    return try std.fmt.allocPrint(b.allocator, "{s}-{s}-{s}", triple);
}
