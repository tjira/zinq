.SHELLFLAGS := $(if $(filter $(OS),Windows_NT),-NoProfile -Command,-c)

CROSS       ?= 0
DEBUG       ?= 0
INCREMENTAL ?= 0
SUMMARY     ?= 0
TIME_REPORT ?= 0
WATCH       ?= 0
WEBUI	    ?= 0

SHELL := $(if $(filter $(OS),Windows_NT),powershell.exe,sh)

ARCH := $(if $(filter $(OS),Windows_NT),x86_64,$(shell uname -m | tr '[:upper:]' '[:lower:]' | sed 's/arm64/aarch64/'))
OS   := $(if $(filter $(OS),Windows_NT),windows,$(shell uname -s | tr '[:upper:]' '[:lower:]' | sed 's/darwin/macos/'))

ZIG_VERSION := 0.15.1
ZLS_VERSION := 0.15.0

ZIG_FLAGS := $(if $(filter 1,$(DEBUG)),-Doptimize=Debug,-Doptimize=ReleaseFast)
ZIG_FLAGS := $(if $(filter 1,$(CROSS)),$(ZIG_FLAGS),$(ZIG_FLAGS) -Dtarget=$(ARCH)-$(OS))
ZIG_FLAGS := $(if $(filter 1,$(INCREMENTAL)),$(ZIG_FLAGS) -fincremental,$(ZIG_FLAGS))
ZIG_FLAGS := $(if $(filter 1,$(SUMMARY)),$(ZIG_FLAGS) --summary all,$(ZIG_FLAGS))
ZIG_FLAGS := $(if $(filter 1,$(TIME_REPORT)),$(ZIG_FLAGS) --time-report,$(ZIG_FLAGS))
ZIG_FLAGS := $(if $(filter 1,$(WATCH)),$(ZIG_FLAGS) --watch,$(ZIG_FLAGS))
ZIG_FLAGS := $(if $(filter 1,$(WEBUI)),$(ZIG_FLAGS) --webui=[::1]:12345,$(ZIG_FLAGS))

all: zinq

# ACORN BUILDING TARGETS ===============================================================================================================================================================================

.PHONY: zinq run test

zinq: .zig-bin/zig$(if $(filter $(OS),windows),.exe)
	./.zig-bin/zig build $(ZIG_FLAGS)

run: .zig-bin/zig$(if $(filter $(OS),windows),.exe)
	./.zig-bin/zig build $(ZIG_FLAGS) run

test: .zig-bin/zig$(if $(filter $(OS),windows),.exe)
	./.zig-bin/zig build $(ZIG_FLAGS) test

# CUSTOM TARGETS =======================================================================================================================================================================================

docs: .zig-bin/zig$(if $(filter $(OS),windows),.exe)
	@./.zig-bin/zig build-lib -femit-docs=docs/code -fno-emit-bin src/main.zig

linguist:
	@github-linguist

profile: zinq
	@valgrind --callgrind-out-file=callgrind.out --tool=callgrind ./zig-out/$(ARCH)-$(OS)/zinq
	@gprof2dot -e 1 -f callgrind -n 5 -o profile.dot -z main.main callgrind.out
	@dot -T pdf profile.dot -o profile.pdf

serve: docs
	@cd docs && bundle exec jekyll serve

wheel:
	python -m build

# EXTERNAL TARGETS =====================================================================================================================================================================================

ifeq ($(OS), windows)
.zig-bin/zig.exe:
	cmd /c curl -Ls -o zig.zip https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).zip
	cmd /c curl -Ls -o zls.zip https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-$(ARCH)-$(OS).zip
	Add-Type -AssemblyName System.IO.Compression.FileSystem ; [System.IO.Compression.ZipFile]::ExtractToDirectory("$$PWD\zig.zip", "$$PWD")
	Move-Item zig-$(ARCH)-$(OS)-$(ZIG_VERSION) .zig-bin ; Remove-Item .zig-bin/LICENSE, .zig-bin/README.md
	Add-Type -AssemblyName System.IO.Compression.FileSystem ; [System.IO.Compression.ZipFile]::ExtractToDirectory("$$PWD\zls.zip", "$$PWD/.zig-bin")
	Remove-Item zig.zip, zls.zip
else
.zig-bin/zig:
	mkdir -p .zig-bin && wget -q -O - https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).tar.xz | tar -Jx -C .zig-bin --strip-components=1 && touch .zig-bin/zig
	mkdir -p .zig-bin && wget -q -O - https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-$(ARCH)-$(OS).tar.xz | tar -Jx -C .zig-bin --strip-components=0 && touch .zig-bin/zls
endif

# CLEAN TARGETS ========================================================================================================================================================================================

clean:
	git clean -dffx
