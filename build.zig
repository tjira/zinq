const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});

    const exe = b.addExecutable(.{
        .name = "zinq",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = b.standardTargetOptions(.{}),
            .strip = optimize != .Debug,
            .optimize = optimize,
            .link_libc = true,
        }),
    });

    const libinc = .{ "lib", "inclde" };

    exe.root_module.addLibraryPath(.{ .cwd_relative = "external-x86_64-linux/" ++ libinc[0] });
    exe.root_module.addIncludePath(.{ .cwd_relative = "external-x86_64-linux/" ++ libinc[1] });

    exe.root_module.linkSystemLibrary("fftw3", .{});

    b.installArtifact(exe);

    const run_exe = b.addRunArtifact(exe);

    b.step("run", "Run the application").dependOn(&run_exe.step);
}
