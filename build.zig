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

    const use_libint   = builder.option(bool, "use-libint",   "Link libint (library for computing molecular integrals)"               ) orelse false;
    const use_openblas = builder.option(bool, "use-openblas", "Link openblas (library with optimized BLAS and LAPACK implementations)") orelse false;
    const use_xc       = builder.option(bool, "use-xc",       "Link libxc (library with various XC functionals)"                      ) orelse false;

    const link_c = use_xc or use_openblas or use_libint;

    const main_options = try generateOptions(builder, use_libint, use_openblas, use_xc);

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

    if (link_c) {

        try addPaths(builder, main_executable, target);

        main_executable.linkLibC(); if (use_libint) main_executable.linkLibCpp();

        if (use_libint) main_executable.root_module.addCSourceFile(.{.file = builder.path("src/libint.cpp"), .flags = &[_][]const u8{}, .language = .cpp});

        if (use_libint  ) main_executable.linkSystemLibrary("int2"    );
        if (use_openblas) main_executable.linkSystemLibrary("openblas");
        if (use_xc      ) main_executable.linkSystemLibrary("xc"      );
    }

    const test_executable = builder.addTest(.{
        .name = "test", .root_module = main_executable.root_module
    });

    const docs_target = builder.addInstallDirectory(.{
        .source_dir = builder.addLibrary(.{
            .name = "main", .root_module = main_executable.root_module,
        }).getEmittedDocs(), .install_dir = .prefix, .install_subdir = "../docs/code"
    });

    main_executable.root_module.addOptions("config", main_options);

    try install(builder, main_executable, builder.getInstallStep(), link_c); try addTools(builder, main_executable, builder.getInstallStep(), link_c);

    builder.step("docs", "Generate documentation"     ).dependOn(&docs_target                            .step);
    builder.step("run",  "Run the compiled executable").dependOn(&builder.addRunArtifact(main_executable).step);
    builder.step("test", "Run unit tests"             ).dependOn(&builder.addRunArtifact(test_executable).step);

    const cross = builder.step("cross", "Cross-compile for all targets");

    for (0..targets.len) |i| {

        const resolved = builder.resolveTargetQuery(targets[i]);

        const matrix_executable = builder.addExecutable(.{
            .name = "zinq",
            .root_module = builder.createModule(.{
                .optimize = optimize,
                .root_source_file = builder.path("src/main.zig"),
                .strip = optimize != .Debug,
                .single_threaded = false,
                .target = resolved
            })
        });

        matrix_executable.root_module.addOptions("config", try generateOptions(builder, false, false, false));

        try install(builder, matrix_executable, cross, false); try addTools(builder, matrix_executable, cross, false);

        if (link_c and checkLibcLinkAvailability(resolved)) try addLibcLinkedExecutables(builder, optimize, targets[i], cross, use_libint, use_openblas, use_xc);
    }
}

pub fn checkLibcLinkAvailability(target: std.Build.ResolvedTarget) bool {
    const os = target.result.os.tag; const arch = target.result.cpu.arch;

    if (os == .linux and arch == .x86_64 ) return true;
    if (os == .linux and arch == .aarch64) return true;

    return false;
}

pub fn addLibcLinkedExecutables(builder: *std.Build, optimize: std.builtin.OptimizeMode, target_query: std.Target.Query, step: *std.Build.Step, use_libint: bool, use_openblas: bool, use_xc: bool) !void {
    const target = builder.resolveTargetQuery(target_query);

    const libc_linked_executable = builder.addExecutable(.{
        .name = "zinq",
        .root_module = builder.createModule(.{
            .optimize = optimize,
            .root_source_file = builder.path("src/main.zig"),
            .strip = optimize != .Debug,
            .single_threaded = false,
            .target = target
        })
    });

    try addPaths(builder, libc_linked_executable, target);

    libc_linked_executable.linkLibC(); if (use_libint) libc_linked_executable.linkLibCpp();

    if (use_libint) libc_linked_executable.root_module.addCSourceFile(.{.file = builder.path("src/libint.cpp"), .flags = &[_][]const u8{}, .language = .cpp});

    if (use_libint  ) libc_linked_executable.linkSystemLibrary("int2"    );
    if (use_openblas) libc_linked_executable.linkSystemLibrary("openblas");
    if (use_xc      ) libc_linked_executable.linkSystemLibrary("xc"      );

    libc_linked_executable.root_module.addOptions("config", try generateOptions(builder, use_libint, use_openblas, use_xc));

    try install(builder, libc_linked_executable, step, true); try addTools(builder, libc_linked_executable, step, true);
}

pub fn addPaths(builder: *std.Build, executable: *std.Build.Step.Compile, target: std.Build.ResolvedTarget) !void {
    const arch_name = @tagName(target.result.cpu.arch);
    const os_name   = @tagName(target.result.os.tag  );

    const include        = try std.mem.concat(builder.allocator, u8, &[_][]const u8{"external-", arch_name, "-", os_name, "/include"        });
    const include_eigen  = try std.mem.concat(builder.allocator, u8, &[_][]const u8{"external-", arch_name, "-", os_name, "/include/eigen3" });
    const include_libint = try std.mem.concat(builder.allocator, u8, &[_][]const u8{"external-", arch_name, "-", os_name, "/include/libint2"});
    const lib            = try std.mem.concat(builder.allocator, u8, &[_][]const u8{"external-", arch_name, "-", os_name, "/lib"            });

    executable.root_module.addIncludePath(.{.cwd_relative = "include"     });
    executable.root_module.addIncludePath(.{.cwd_relative = include       });
    executable.root_module.addIncludePath(.{.cwd_relative = include_eigen });
    executable.root_module.addIncludePath(.{.cwd_relative = include_libint});
    executable.root_module.addLibraryPath(.{.cwd_relative = lib           });
}

pub fn addTools(builder: *std.Build, main_executable: *std.Build.Step.Compile, step: *std.Build.Step, libc: bool) !void {
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

        try install(builder, tool_executable, step, libc);
    }
}

pub fn generateOptions(builder: *std.Build, use_libint: bool, use_openblas: bool, use_xc: bool) !*std.Build.Step.Options {
    const options = builder.addOptions();

    options.addOption([]const u8, "libint_version",   try getVersion(builder, "lib/libint"  ));
    options.addOption([]const u8, "openblas_version", try getVersion(builder, "lib/openblas"));
    options.addOption([]const u8, "libxc_version",    try getVersion(builder, "lib/libxc"   ));
    options.addOption([]const u8, "zinq_version",     try getVersion(builder, null          ));

    options.addOption(bool, "use_libint",   use_libint  );
    options.addOption(bool, "use_openblas", use_openblas);
    options.addOption(bool, "use_xc",       use_xc      );

    return options;
}

fn getVersion(builder: *std.Build, dir: ?[]const u8) ![]const u8 {
    const cwd = if (dir) |d| std.fs.cwd().openDir(d, .{}) catch return "UNKNOWN" else std.fs.cwd();

    const result = std.process.Child.run(.{.allocator = builder.allocator, .argv = &.{"git", "describe", "--tags"}, .cwd_dir = cwd}) catch {
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

pub fn install(builder: *std.Build, main_executable: *std.Build.Step.Compile, step: *std.Build.Step, libc: bool) !void {
    const target = main_executable.root_module.resolved_target.?;

    if (target.query.isNative()) {

        const main_executable_install = builder.addInstallArtifact(main_executable, .{});

        step.dependOn(&main_executable_install.step); return;
    }

    const dest = try std.fmt.allocPrint(builder.allocator, "{s}-{s}{s}", .{
        @tagName(target.result.cpu.arch), @tagName(target.result.os.tag), if (libc) "-libc" else ""
    });

    const main_executable_install = builder.addInstallArtifact(main_executable, .{
        .dest_dir = .{.override = .{.custom = dest}}
    });

    step.dependOn(&main_executable_install.step);
}
