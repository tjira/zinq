const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});

    const zinq_module = setupZinq(b, optimize, b.standardTargetOptions(.{}));

    setupTests(b, zinq_module);
}

fn setupZinq(b: *std.Build, opt: std.builtin.OptimizeMode, target: std.Build.ResolvedTarget) *std.Build.Module {
    const zinq_module = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = opt,
        .strip = opt == .ReleaseFast or opt == .ReleaseSmall,
        .link_libc = true,
    });

    linkDependencies(zinq_module);

    const exe_zinq = b.addExecutable(.{
        .name = "zinq",
        .root_module = zinq_module,
    });

    b.installArtifact(exe_zinq);

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
    });

    test_module.addImport("zinq", zinq_module);

    const exe_test = b.addTest(.{
        .root_module = test_module,
    });

    const run_exe_test = b.addRunArtifact(exe_test);

    b.step("test", "Run unit tests").dependOn(&run_exe_test.step);
}

fn linkDependencies(module: *std.Build.Module) void {
    const dirs = [_][]const u8{ "lib", "include" };

    module.addLibraryPath(.{ .cwd_relative = "external-x86_64-linux/" ++ dirs[0] });
    module.addIncludePath(.{ .cwd_relative = "external-x86_64-linux/" ++ dirs[1] });

    const libs = [_][]const u8{
        "fftw3",
        "openblas",
    };

    for (libs) |lib| {
        module.linkSystemLibrary(lib, .{});
    }
}
