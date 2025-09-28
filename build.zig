const builtin = @import("builtin");
const std = @import("std");

const targets: []const std.Target.Query = &.{
    .{.os_tag = .freebsd, .cpu_arch = .aarch64},
    .{.os_tag = .freebsd, .cpu_arch = .x86_64 },
    .{.os_tag = .linux,   .cpu_arch = .aarch64},
    .{.os_tag = .linux,   .cpu_arch = .arm    },
    .{.os_tag = .linux,   .cpu_arch = .riscv32},
    .{.os_tag = .linux,   .cpu_arch = .riscv64},
    .{.os_tag = .linux,   .cpu_arch = .x86    },
    .{.os_tag = .linux,   .cpu_arch = .x86_64 },
    .{.os_tag = .macos,   .cpu_arch = .aarch64},
    .{.os_tag = .macos,   .cpu_arch = .x86_64 },
    .{.os_tag = .netbsd,  .cpu_arch = .aarch64},
    .{.os_tag = .netbsd,  .cpu_arch = .x86_64 },
    .{.os_tag = .windows, .cpu_arch = .aarch64},
    .{.os_tag = .windows, .cpu_arch = .x86_64 }
};

pub fn build(builder: *std.Build) !void {
    const optimize = builder.standardOptimizeOption(.{});
    const target = builder.standardTargetOptions(.{});

    for (0..targets.len) |i| {

        const target_query = if (target.query.cpu_arch == null) targets[i] else target.query;

        const main_executable = builder.addExecutable(.{
            .name = "zinq",
            .root_module = builder.createModule(.{
                .optimize = optimize,
                .root_source_file = builder.path("src/main.zig"),
                .strip = optimize != .Debug,
                .single_threaded = false,
                .target = builder.resolveTargetQuery(target_query)
            }),
            .use_llvm = true // needed for valgrind for now
        });

        const main_executable_install = builder.addInstallArtifact(main_executable, .{
            .dest_dir = .{.override = .{.custom = try target_query.zigTriple(builder.allocator)}}
        });

        builder.getInstallStep().dependOn(&main_executable_install.step);

        if (builtin.target.cpu.arch == target_query.cpu_arch and builtin.target.os.tag == target_query.os_tag) {

            const test_executable = builder.addTest(.{
                .name = "test", .root_module = main_executable.root_module,
            });

            builder.step("run",  "Run the compiled executable").dependOn(&builder.addRunArtifact(main_executable).step);
            builder.step("test", "Run unit tests"             ).dependOn(&builder.addRunArtifact(test_executable).step);
        }

        if (target.query.cpu_arch != null) break;
    }
}
