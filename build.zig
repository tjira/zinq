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

    const use_xc = builder.option(bool, "use-xc", "Link libxc (library with various XC functionals)") orelse false;

    const options = generateOptions(builder, use_xc);

    const main_executable = builder.addExecutable(.{
        .name = "zinq",
        .root_module = builder.createModule(.{
            .optimize = optimize,
            .root_source_file = builder.path("src/main.zig"),
            .strip = optimize != .Debug,
            .single_threaded = false,
            .target = target
        })
    });

    if (use_xc) {
        main_executable.root_module.addIncludePath(.{.cwd_relative = "external/include"});
        main_executable.root_module.addLibraryPath(.{.cwd_relative = "external/lib"    });
    }

    if (use_xc) {
        main_executable.linkLibC(); main_executable.linkSystemLibrary("xc");
    }

    const test_executable = builder.addTest(.{
        .name = "test", .root_module = main_executable.root_module
    });

    const docs_target = builder.addInstallDirectory(.{
        .source_dir = builder.addLibrary(.{
            .name = "main", .root_module = main_executable.root_module,
        }).getEmittedDocs(), .install_dir = .prefix, .install_subdir = "../docs/code"
    });

    main_executable.root_module.addOptions("config", options);

    try install(builder, main_executable, builder.getInstallStep()); try addTools(builder, main_executable, builder.getInstallStep());

    builder.step("docs", "Generate documentation"     ).dependOn(&docs_target                            .step);
    builder.step("run",  "Run the compiled executable").dependOn(&builder.addRunArtifact(main_executable).step);
    builder.step("test", "Run unit tests"             ).dependOn(&builder.addRunArtifact(test_executable).step);

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
            })
        });

        matrix_executable.root_module.addOptions("config", options);

        try install(builder, matrix_executable, cross); try addTools(builder, matrix_executable, cross);
    }
}

pub fn addTools(builder: *std.Build, main_executable: *std.Build.Step.Compile, step: *std.Build.Step) !void {
    const optimize = main_executable.root_module.optimize;
    const target = main_executable.root_module.resolved_target.?;

    var tool_dir = try std.fs.cwd().openDir("tool", .{.iterate = true}); defer tool_dir.close();

    var iterator = tool_dir.iterate();

    while (try iterator.next()) |entry| {

        const tool = entry.name[0..entry.name.len - 4];

        const tool_executable = builder.addExecutable(.{
            .name = builder.fmt("zinq-{s}", .{tool}),
            .root_module = builder.createModule(.{
                .optimize = optimize,
                .root_source_file = builder.path(builder.fmt("tool/{s}.zig", .{tool})),
                .strip = optimize != .Debug,
                .single_threaded = true,
                .target = target
            })
        });

        tool_executable.root_module.addImport("zinq", main_executable.root_module);

        try install(builder, tool_executable, step);
    }
}

pub fn generateOptions(builder: *std.Build, use_xc: bool) *std.Build.Step.Options {
    const options = builder.addOptions();

    options.addOption([]const u8, "zinq_version", getVersion(builder));

    options.addOption(bool, "use_xc", use_xc);

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

pub fn install(builder: *std.Build, main_executable: *std.Build.Step.Compile, step: *std.Build.Step) !void {
    const target = main_executable.root_module.resolved_target.?;

    if (target.query.isNative()) {

        const main_executable_install = builder.addInstallArtifact(main_executable, .{});

        step.dependOn(&main_executable_install.step); return;
    }

    const dest = try std.fmt.allocPrint(builder.allocator, "{s}-{s}", .{
        @tagName(target.result.cpu.arch), @tagName(target.result.os.tag)
    });

    const main_executable_install = builder.addInstallArtifact(main_executable, .{
        .dest_dir = .{.override = .{.custom = dest}}
    });

    step.dependOn(&main_executable_install.step);
}
