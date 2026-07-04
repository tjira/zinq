DEBUG ?= 0

ARCH := $(if $(filter $(OS),Windows_NT),x86_64,$(shell uname -m | tr '[:upper:]' '[:lower:]' | sed 's/arm64/aarch64/'))
OS   := $(if $(filter $(OS),Windows_NT),windows,$(shell uname -s | tr '[:upper:]' '[:lower:]' | sed 's/darwin/macos/'))

ZIG_VERSION := 0.16.0
ZLS_VERSION := 0.16.0

.PHONY: all zinq run test

all: zinq

zinq:
	@zig build $(if $(filter 0,$(DEBUG)),--release=fast)

run:
	@zig build $(if $(filter 0,$(DEBUG)),--release=fast) run

test:
	@zig build $(if $(filter 0,$(DEBUG)),--release=fast) test

library:
	@curl -Ls https://nightly.link/tjira/zinq/workflows/library/master/external-$(ARCH)-$(OS)-musl.zip | bsdtar -xf -

zig:
	@mkdir -p .zig-bin && curl -Ls https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).tar.xz | tar -Jx -C .zig-bin --strip-components=1

zls:
	@mkdir -p .zig-bin && curl -Ls https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-$(ARCH)-$(OS).tar.xz | tar -Jx -C .zig-bin --strip-components=0

clean:
	@git clean -dffx
