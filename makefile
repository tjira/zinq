.SHELLFLAGS := $(if $(filter $(OS),Windows_NT),-NoProfile -Command)

DEBUG   ?= 0
SUMMARY ?= 0

ARCH := $(if $(filter $(OS),Windows_NT),x86_64,$(shell uname -m | tr '[:upper:]' '[:lower:]' | sed 's/arm64/aarch64/'))
OS   := $(if $(filter $(OS),Windows_NT),windows,$(shell uname -s | tr '[:upper:]' '[:lower:]' | sed 's/darwin/macos/'))

SHELL := $(if $(filter $(OS),windows),powershell.exe,bash)

ZIG_VERSION := 0.15.1
ZLS_VERSION := 0.15.0

ZIG_FLAGS := --summary $(if $(filter 1,$(SUMMARY)),all,none) $(if $(filter 1,$(DEBUG)),-Doptimize=Debug,-Doptimize=ReleaseFast)

all: zinq

# ACORN BUILDING TARGETS ===============================================================================================================================================================================

.PHONY: zinq run test

zinq: .zig-bin/zig$(if $(filter $(OS),windows),.exe) clean-output
	./.zig-bin/zig build $(ZIG_FLAGS)

run: .zig-bin/zig$(if $(filter $(OS),windows),.exe)
	./.zig-bin/zig build $(ZIG_FLAGS) run

test: .zig-bin/zig$(if $(filter $(OS),windows),.exe)
	./.zig-bin/zig build $(ZIG_FLAGS) test

# CUSTOM TARGETS =======================================================================================================================================================================================

docs: .zig-bin/zig$(if $(filter $(OS),windows),.exe) clean-docs
	@./.zig-bin/zig build-lib -femit-docs=docs/code -fno-emit-bin src/main.zig

linguist:
	@github-linguist

profile: zinq
	@valgrind --callgrind-out-file=callgrind.out --tool=callgrind ./zig-out/$(ARCH)-$(OS)/zinq
	@gprof2dot -e 1 -f callgrind -n 5 -o profile.dot -z main.main callgrind.out
	@dot -T pdf profile.dot -o profile.pdf

serve: docs
	@cd docs && bundle exec jekyll serve

# EXTERNAL TARGETS =====================================================================================================================================================================================

ifeq ($(OS), windows)
.zig-bin/zig.exe:
	curl -o zig.zip https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).zip
	curl -o zls.zip https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-$(ARCH)-$(OS).zip
	Expand-Archive -Force -Path zig.zip -DestinationPath . ; mv zig-$(ARCH)-$(OS)-$(ZIG_VERSION) .zig-bin ; Expand-Archive -Force -Path zls.zip -DestinationPath .zig-bin
	Remove-Item zig.zip, zls.zip
else
.zig-bin/zig:
	mkdir -p .zig-bin && wget -q -O - https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).tar.xz | tar -Jx -C .zig-bin --strip-components=1 && touch .zig-bin/zig
	mkdir -p .zig-bin && wget -q -O - https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-$(ARCH)-$(OS).tar.xz | tar -Jx -C .zig-bin --strip-components=0 && touch .zig-bin/zls
endif

# CLEAN TARGETS ========================================================================================================================================================================================

clean: clean-archive clean-docs clean-output clean-root clean-zig

ifeq ($(OS), windows)
clean-archive:
	@&{Remove-Item -ErrorAction SilentlyContinue -Force -Recurse archive/*.mat, archive/*.ppm}
else
clean-archive:
	@rm -rf archive/*.mat archive/*.ppm
endif

ifeq ($(OS), windows)
clean-docs:
	@&{Remove-Item -ErrorAction SilentlyContinue -Force -Recurse docs/_site, docs/.jekyll-cache, docs/code, docs/tex, docs/*.locked}
else
clean-docs:
	@rm -rf docs/_site docs/.jekyll-cache docs/code docs/tex docs/*.locked
endif

ifeq ($(OS), windows)
clean-output:
	@&{Remove-Item -ErrorAction SilentlyContinue -Force -Recurse zig-out}
else
clean-output:
	@rm -rf zig-out
endif

ifeq ($(OS), windows)
clean-root:
	@&{Remove-Item -ErrorAction SilentlyContinue -Force -Recurse *.all, *.dot, *.gz, *.json, *.mat, *.out, *.pdf, *.png, *.ppm, *.xyz, *.zip}
else
clean-root:
	@rm -rf *.all *.dot *.gz *.json *.mat *.out *.pdf *.png *.ppm *.xyz *.zip
endif

ifeq ($(OS), windows)
clean-zig:
	@&{Remove-Item -ErrorAction SilentlyContinue -Force -Recurse ${HOME}/.cache/zig, .zig-bin, .zig-cache}
else
clean-zig:
	@rm -rf ${HOME}/.cache/zig .zig-bin .zig-cache
endif
