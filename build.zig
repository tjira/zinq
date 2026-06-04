const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});

    const zinq_module = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .link_libc = true,
        .optimize = optimize,
        .strip = optimize != .Debug,
        .target = b.standardTargetOptions(.{}),
    });

    link(zinq_module);

    const exe = b.addExecutable(.{
        .name = "zinq",
        .root_module = zinq_module,
    });

    b.installArtifact(exe);

    const run_exe = b.addRunArtifact(exe);

    b.step("run", "Run the application").dependOn(&run_exe.step);

    const exe_test = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("test/main.zig"),
            .target = exe.root_module.resolved_target,
            .strip = exe.root_module.strip,
            .optimize = exe.root_module.optimize,
            .link_libc = exe.root_module.link_libc,
        }),
    });

    link(exe_test.root_module);

    exe_test.root_module.addImport("zinq", zinq_module);

    const run_exe_test = b.addRunArtifact(exe_test);

    b.step("test", "Run unit tests").dependOn(&run_exe_test.step);
}

pub fn link(module: *std.Build.Module) void {
    // zig fmt: off
    module.addLibraryPath(.{ .cwd_relative = "external-x86_64-linux/lib"     });
    module.addIncludePath(.{ .cwd_relative = "external-x86_64-linux/include" });
    // zig fmt: on

    // zig fmt: off
    module.linkSystemLibrary("fftw3",    .{});
    module.linkSystemLibrary("openblas", .{});
    // zig fmt: on
}
