.SHELLFLAGS := $(if $(filter $(OS),Windows_NT),-NoProfile -Command,-c)

DEBUG ?= 0

SHELL := $(if $(filter $(OS),Windows_NT),powershell.exe,sh)

ARCH := $(if $(filter $(OS),Windows_NT),x86_64,$(shell uname -m | tr '[:upper:]' '[:lower:]' | sed 's/arm64/aarch64/'))
OS   := $(if $(filter $(OS),Windows_NT),windows,$(shell uname -s | tr '[:upper:]' '[:lower:]' | sed 's/darwin/macos/'))

ZIG_VERSION := 0.16.0
ZLS_VERSION := 0.16.0

HAS_ZIG := $(shell $(if $(filter windows,$(OS)),(Get-Command zig -ErrorAction SilentlyContinue).Path,command -v zig 2> /dev/null))
HAS_ZLS := $(shell $(if $(filter windows,$(OS)),(Get-Command zls -ErrorAction SilentlyContinue).Path,command -v zls 2> /dev/null))

COMPILER := $(if $(HAS_ZIG),zig,./.zig-bin/zig$(if $(filter $(OS),windows),.exe))

.PHONY: all zinq docs run test

all: .env.fish .env.ps1 .env.sh zinq

# ZINQ BUILDING TARGETS ========================================================================================================================================

zinq: $(if $(HAS_ZIG),,.zig-bin/zig$(if $(filter $(OS),windows),.exe)) $(if $(HAS_ZLS),,.zig-bin/zls$(if $(filter $(OS),windows),.exe)) external-$(ARCH)-$(OS)
	@$(COMPILER) build $(if $(filter 0,$(DEBUG)),--release=fast)

docs: $(if $(HAS_ZIG),,.zig-bin/zig$(if $(filter $(OS),windows),.exe)) $(if $(HAS_ZLS),,.zig-bin/zls$(if $(filter $(OS),windows),.exe)) external-$(ARCH)-$(OS)
	@$(COMPILER) build docs

run: $(if $(HAS_ZIG),,.zig-bin/zig$(if $(filter $(OS),windows),.exe)) $(if $(HAS_ZLS),,.zig-bin/zls$(if $(filter $(OS),windows),.exe)) external-$(ARCH)-$(OS)
	@$(COMPILER) build $(if $(filter 0,$(DEBUG)),--release=fast) run

test: $(if $(HAS_ZIG),,.zig-bin/zig$(if $(filter $(OS),windows),.exe)) $(if $(HAS_ZLS),,.zig-bin/zls$(if $(filter $(OS),windows),.exe)) external-$(ARCH)-$(OS)
	@$(COMPILER) build $(if $(filter 0,$(DEBUG)),--release=fast) test

# LIBRARY INSTALLATION TARGETS =================================================================================================================================

external-$(ARCH)-$(OS):
	@curl -Ls https://nightly.link/tjira/zinq/workflows/library/master/external-$(ARCH)-$(OS).zip | bsdtar -xf -

# ENVIRONMENT SCRIPTS ==========================================================================================================================================

.env.fish:
	@echo "fish_add_path $(CURDIR)/.zig-bin" > .env.fish

.env.ps1:
	@echo '$$env:PATH="$(CURDIR)/.zig-bin;$$env:PATH"' > .env.ps1

.env.sh:
	@echo "export PATH=$(CURDIR)/.zig-bin:\$$PATH" > .env.sh

# COMPILER AND LANGUAGE SERVER INSTALLATION TARGETS ============================================================================================================

ifeq ($(OS),windows)
.zig-bin/zig.exe: | .zig-bin
	@curl.exe -Ls -o zig.zip https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).zip ; tar -xf zig.zip -C .zig-bin --strip-components=1 ; rm zig.zip
else
.zig-bin/zig: | .zig-bin
	@curl -Ls https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).tar.xz | tar -Jx -C .zig-bin --strip-components=1
endif

ifeq ($(OS),windows)
.zig-bin/zls.exe: | .zig-bin
	@curl.exe -Ls -o zls.zip https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-$(ARCH)-$(OS).zip ; tar -xf zls.zip -C .zig-bin --strip-components=0 ; rm zls.zip
else
.zig-bin/zls: | .zig-bin
	@curl -Ls https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-$(ARCH)-$(OS).tar.xz | tar -Jx -C .zig-bin --strip-components=0
endif

# DIRECTORY CREATION TARGETS ===================================================================================================================================

.zig-bin:
	@$(if $(filter windows,$(OS)),mkdir .zig-bin -Force | Out-Null,mkdir -p .zig-bin)

# ADDITIONAL TARGETS ===========================================================================================================================================

clean:
	@git clean -dffx
