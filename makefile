.RECIPEPREFIX = >

DEBUG   ?= 0
SUMMARY ?= 0

ARCH := $(shell uname -m | tr '[:upper:]' '[:lower:]')
OS   := $(shell uname -s | tr '[:upper:]' '[:lower:]')

ZIG_VERSION := 0.15.1
ZLS_VERSION := 0.15.0

ZIG_FLAGS := --summary $(if $(filter 1,$(SUMMARY)),all,none) $(if $(filter 1,$(DEBUG)),-Doptimize=Debug,-Doptimize=ReleaseFast)

all: zinq

# ACORN BUILDING TARGETS ===============================================================================================================================================================================

.PHONY: acorn run test

zinq: zig-bin/zig clean-output
> ./zig-bin/zig build $(ZIG_FLAGS)

run: zig-bin/zig
> ./zig-bin/zig build $(ZIG_FLAGS) run

test: zig-bin/zig
> ./zig-bin/zig build $(ZIG_FLAGS) test

# CUSTOM TARGETS =======================================================================================================================================================================================

docs: zig-bin/zig clean-docs
> @./zig-bin/zig build-lib -femit-docs=docs/code -fno-emit-bin src/main.zig

linguist:
> @github-linguist

profile: zinq
> @valgrind --callgrind-out-file=callgrind.out --tool=callgrind ./zig-out/$(ARCH)-$(OS)/zinq
> @gprof2dot --format=callgrind --output=profile.dot --root=main.main callgrind.out
> @dot -T pdf profile.dot -o profile.pdf

serve: docs
> @cd docs && bundle exec jekyll serve

# EXTERNAL TARGETS =====================================================================================================================================================================================

zig-bin/zig:
> mkdir -p zig-bin && wget -q -O - https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).tar.xz | tar -Jx -C zig-bin --strip-components=1 && touch zig-bin/zig
> mkdir -p zig-bin && wget -q -O - https://github.com/zigtools/zls/releases/download/$(ZLS_VERSION)/zls-$(ARCH)-$(OS).tar.xz | tar -Jx -C zig-bin --strip-components=0 && touch zig-bin/zls

# CLEAN TARGETS ========================================================================================================================================================================================

clean: clean-cache clean-docs clean-output clean-root clean-zig

clean-cache:
> @rm -rf ${HOME}/.cache/zig .zig-cache

clean-docs:
> @rm -rf docs/_site docs/.jekyll-cache docs/code docs/tex docs/*.locked

clean-output:
> @rm -rf zig-out

clean-root:
> @rm -rf *.all *.dot *.json *.mat *.out *.pdf *.png *.xyz

clean-zig:
> @rm -rf zig-bin
