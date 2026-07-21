DEBUG ?= 0

ARCH := $(if $(filter $(OS),Windows_NT),x86_64,$(shell uname -m | tr '[:upper:]' '[:lower:]' | sed 's/arm64/aarch64/'))
OS   := $(if $(filter $(OS),Windows_NT),windows,$(shell uname -s | tr '[:upper:]' '[:lower:]' | sed 's/darwin/macos/'))

ZIG_VERSION := 0.16.0
ZLS_VERSION := 0.16.0

HAS_ZIG := $(shell command -v zig 2> /dev/null)
HAS_ZLS := $(shell command -v zls 2> /dev/null)
COMPILER := $(if $(HAS_ZIG),zig,./.zig-bin/zig)

.PHONY: all zinq docs run test

all: .env.fish .env.sh zinq

# ZINQ BUILDING TARGETS ========================================================================================================================================

zinq: $(if $(HAS_ZIG),,.zig-bin/zig) $(if $(HAS_ZLS),,.zig-bin/zls) external-$(ARCH)-$(OS)
	@$(COMPILER) build $(if $(filter 0,$(DEBUG)),--release=fast)

docs: $(if $(HAS_ZIG),,.zig-bin/zig) $(if $(HAS_ZLS),,.zig-bin/zls) external-$(ARCH)-$(OS)
	@$(COMPILER) build docs

run: $(if $(HAS_ZIG),,.zig-bin/zig) $(if $(HAS_ZLS),,.zig-bin/zls) external-$(ARCH)-$(OS)
	@$(COMPILER) build $(if $(filter 0,$(DEBUG)),--release=fast) run

test: $(if $(HAS_ZIG),,.zig-bin/zig) $(if $(HAS_ZLS),,.zig-bin/zls) external-$(ARCH)-$(OS)
	@$(COMPILER) build $(if $(filter 0,$(DEBUG)),--release=fast) -Dtarget=native-native-musl test

# FILE OR DIRECTORY TARGETS ====================================================================================================================================

external-$(ARCH)-$(OS):
	@curl -Ls https://nightly.link/tjira/zinq/workflows/library/master/external-$(ARCH)-$(OS).zip | bsdtar -xf -

.env.fish:
	@echo "fish_add_path $(CURDIR)/.zig-bin" > .env.fish

.env.sh:
	@echo "export PATH=$(CURDIR)/.zig-bin:\$$PATH" > .env.sh

.zig-bin/zig:
	@mkdir -p .zig-bin && curl -Ls https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).tar.xz | tar -Jx -C .zig-bin --strip-components=1 && touch $@

.zig-bin/zls:
	@mkdir -p .zig-bin && curl -Ls https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-$(ARCH)-$(OS).tar.xz | tar -Jx -C .zig-bin --strip-components=0 && touch $@

# ADDITIONAL TARGETS ===========================================================================================================================================

clean:
	@git clean -dffx
