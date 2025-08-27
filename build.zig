const builtin = @import("builtin");
const std     = @import("std");

const targets: []const std.Target.Query = &.{
    .{.os_tag = .linux,   .cpu_arch = .aarch64},
    .{.os_tag = .linux,   .cpu_arch = .riscv64},
    .{.os_tag = .linux,   .cpu_arch = .x86_64 },
    .{.os_tag = .macos,   .cpu_arch = .aarch64},
    .{.os_tag = .macos,   .cpu_arch = .x86_64 },
    .{.os_tag = .windows, .cpu_arch = .aarch64},
    .{.os_tag = .windows, .cpu_arch = .x86_64 },
};

pub fn build(builder: *std.Build) !void {
    const debug = builder.option(bool, "DEBUG", "Build everything in the debug mode") orelse false;

    for (targets) |target| {

        const main_executable = builder.addExecutable(.{
            .name = "zinq",
            .root_module = builder.createModule(.{
                .optimize = if (debug) .Debug else .ReleaseFast,
                .root_source_file = builder.path("src/main.zig"),
                .strip = !debug,
                .target = builder.resolveTargetQuery(target)
            })
        });

        const test_executable = builder.addTest(.{
            .name = "test",
            .root_module = builder.createModule(.{
                .optimize = if (debug) .Debug else .ReleaseFast,
                .root_source_file = builder.path("src/main.zig"),
                .strip = !debug,
                .target = builder.resolveTargetQuery(target)
            })
        });

        const main_executable_install = builder.addInstallArtifact(main_executable, .{
            .dest_dir = .{.override = .{.custom = try target.zigTriple(builder.allocator)}}
        });

        builder.getInstallStep().dependOn(&main_executable_install.step);

        if (builtin.target.cpu.arch == target.cpu_arch and builtin.target.os.tag == target.os_tag) {

            builder.step("run",  "Run the compiled executable").dependOn(&builder.addRunArtifact(main_executable).step);
            builder.step("test", "Run unit tests"             ).dependOn(&builder.addRunArtifact(test_executable).step);

            const docs = builder.addInstallDirectory(.{
                .install_dir = .{
                    .custom = "../docs"
                },
                .install_subdir = "code",
                .source_dir = main_executable.getEmittedDocs()
            });

            builder.step("docs", "Install the code documentation").dependOn(&docs.step);
        }
    }
}
