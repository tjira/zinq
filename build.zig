const builtin = @import("builtin");
const std = @import("std");

const targets: []const std.Target.Query = &.{
    .{.os_tag = .freebsd, .cpu_arch = .arm    },
    .{.os_tag = .freebsd, .cpu_arch = .aarch64},
    .{.os_tag = .freebsd, .cpu_arch = .x86    },
    .{.os_tag = .freebsd, .cpu_arch = .x86_64 },
    .{.os_tag = .linux,   .cpu_arch = .arm    },
    .{.os_tag = .linux,   .cpu_arch = .aarch64},
    .{.os_tag = .linux,   .cpu_arch = .riscv64},
    .{.os_tag = .linux,   .cpu_arch = .x86    },
    .{.os_tag = .linux,   .cpu_arch = .x86_64 },
    .{.os_tag = .macos,   .cpu_arch = .aarch64},
    .{.os_tag = .macos,   .cpu_arch = .x86_64 },
    .{.os_tag = .netbsd,  .cpu_arch = .arm    },
    .{.os_tag = .netbsd,  .cpu_arch = .aarch64},
    .{.os_tag = .netbsd,  .cpu_arch = .x86    },
    .{.os_tag = .netbsd,  .cpu_arch = .x86_64 },
    .{.os_tag = .windows, .cpu_arch = .aarch64},
    .{.os_tag = .windows, .cpu_arch = .x86    },
    .{.os_tag = .windows, .cpu_arch = .x86_64 }
};

pub fn build(builder: *std.Build) !void {
    const optimize = builder.standardOptimizeOption(.{});
    const target = builder.standardTargetOptions(.{});

    const options = generateOptions(builder);

    const main_executable = builder.addExecutable(.{
        .name = "zinq",
        .root_module = builder.createModule(.{
            .optimize = optimize,
            .root_source_file = builder.path("src/main.zig"),
            .strip = optimize != .Debug,
            .single_threaded = false,
            .target = target
        }),
        .use_llvm = true // needed for valgrind for now
    });

    const benchmark_executable = builder.addExecutable(.{
        .name = "benchmark",
        .root_module = builder.createModule(.{
            .optimize = optimize,
            .root_source_file = builder.path("tool/benchmark.zig"),
            .strip = optimize != .Debug,
            .single_threaded = true,
            .target = target
        }),
        .use_llvm = true // needed for valgrind for now
    });

    const test_executable = builder.addTest(.{
        .name = "test", .root_module = main_executable.root_module,
    });

    const docs_target = builder.addInstallDirectory(.{
        .source_dir = builder.addLibrary(.{
            .name = "main", .root_module = main_executable.root_module,
        }).getEmittedDocs(), .install_dir = .prefix, .install_subdir = "../docs/code"
    });

    main_executable.root_module.addOptions("config", options);

    benchmark_executable.root_module.addImport("zinq", main_executable.root_module);

    const main_executable_install = builder.addInstallArtifact(main_executable, .{});

    builder.getInstallStep().dependOn(&main_executable_install.step);

    builder.step("benchmark", "Run the compiled executable").dependOn(&builder.addRunArtifact(benchmark_executable).step);
    builder.step("docs",      "Generate documentation"     ).dependOn(&docs_target                                 .step);
    builder.step("run",       "Run the compiled executable").dependOn(&builder.addRunArtifact(main_executable     ).step);
    builder.step("test",      "Run unit tests"             ).dependOn(&builder.addRunArtifact(test_executable     ).step);

    const cross = builder.step("cross", "Cross-compile for all targets");

    for (0..targets.len) |i| {

        const matrix_executable = builder.addExecutable(.{
            .name = "zinq",
            .root_module = builder.createModule(.{
                .optimize = optimize,
                .root_source_file = builder.path("src/main.zig"),
                .strip = optimize != .Debug,
                .single_threaded = false,
                .target = builder.resolveTargetQuery(targets[i])
            }),
            .use_llvm = true // needed for valgrind for now
        });

        matrix_executable.root_module.addOptions("config", options);

        const matrix_executable_install = builder.addInstallArtifact(matrix_executable, .{
            .dest_dir = .{.override = .{.custom = try targets[i].zigTriple(builder.allocator)}}
        });

        cross.dependOn(&matrix_executable_install.step);
    }
}

pub fn generateOptions(builder: *std.Build) *std.Build.Step.Options {
    const options = builder.addOptions();

    options.addOption([]const u8, "zinq_version", getVersion(builder));

    return options;
}

fn getVersion(builder: *std.Build) []const u8 {
    const result = std.process.Child.run(.{.allocator = builder.allocator, .argv = &.{"git", "describe", "--tags"}}) catch {
        return "UNKNOWN";
    };

    defer {
        builder.allocator.free(result.stdout);
        builder.allocator.free(result.stderr);
    }

    if (result.term.Exited != 0) {
        return "UNKNOWN";
    }

    const version = builder.allocator.dupe(u8, std.mem.trim(u8, result.stdout, " \n\r\t")) catch return "UNKNOWN";

    for (version) |*c| if (c.* == '-') {
        c.* = '+'; break;
    };

    for (version) |*c| if (c.* == '-') {
        c.* = '.'; break;
    };

    return version;
}
