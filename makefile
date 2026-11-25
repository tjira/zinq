.SHELLFLAGS := $(if $(filter $(OS),Windows_NT),-NoProfile -Command,-c)

ARGS  ?=
NODE  ?=

DEBUG   ?= 0
GENERIC ?= 0

SHELL := $(if $(filter $(OS),Windows_NT),powershell.exe,sh)

ARCH := $(if $(filter $(OS),Windows_NT),x86_64,$(shell uname -m | tr '[:upper:]' '[:lower:]' | sed 's/arm64/aarch64/'))
OS   := $(if $(filter $(OS),Windows_NT),windows,$(shell uname -s | tr '[:upper:]' '[:lower:]' | sed 's/darwin/macos/'))

ZIG_VERSION := 0.17.1

HAS_ZIG := $(shell $(if $(filter windows,$(OS)),(Get-Command zig -ErrorAction SilentlyContinue).Path,command -v zig 2> /dev/null))

COMPILER := $(if $(HAS_ZIG),zig,./.zig-bin/zig$(if $(filter $(OS),windows),.exe))

.PHONY: all basis zinq docs fmt molecule run test $(if $(filter linux,$(OS)),callgrind profile)

all: .env.fish .env.ps1 .env.sh zinq

# ZINQ BUILDING TARGETS ========================================================================================================================================

zinq: $(if $(HAS_ZIG),,.zig-bin/zig$(if $(filter $(OS),windows),.exe)) external-$(ARCH)-$(OS)
	@$(COMPILER) build $(if $(filter 0,$(DEBUG)),--release=fast)

docs: $(if $(HAS_ZIG),,.zig-bin/zig$(if $(filter $(OS),windows),.exe)) external-$(ARCH)-$(OS)
	@$(COMPILER) build docs

fmt: $(if $(HAS_ZIG),,.zig-bin/zig$(if $(filter $(OS),windows),.exe)) external-$(ARCH)-$(OS)
	@$(COMPILER) fmt src test build.zig

run: $(if $(HAS_ZIG),,.zig-bin/zig$(if $(filter $(OS),windows),.exe)) external-$(ARCH)-$(OS)
	@$(COMPILER) build $(if $(filter 0,$(DEBUG)),--release=fast) run $(if $(strip $(ARGS)),-- $(ARGS))

test: $(if $(HAS_ZIG),,.zig-bin/zig$(if $(filter $(OS),windows),.exe)) external-$(ARCH)-$(OS)
	@$(COMPILER) build $(if $(filter 0,$(DEBUG)),--release=fast) test

# LIBRARY INSTALLATION TARGETS =================================================================================================================================

ifeq ($(OS),windows)
external-$(ARCH)-$(OS):
	@curl.exe -L# -o external.zip https://nightly.link/tjira/zinq/workflows/library/master/external-$(ARCH)-$(OS)$(if $(filter 1,$(GENERIC)),-generic).zip ; tar -xf external.zip ; rm external.zip
else
external-$(ARCH)-$(OS):
	@curl -L# -o external.zip https://nightly.link/tjira/zinq/workflows/library/master/external-$(ARCH)-$(OS)$(if $(filter 1,$(GENERIC)),-generic).zip ; unzip -q external.zip ; rm external.zip
endif

# ENVIRONMENT SCRIPTS ==========================================================================================================================================

.env.fish:
	@echo "fish_add_path $(CURDIR)/.zig-bin" > .env.fish

.env.ps1:
	@echo '$$env:PATH="$(CURDIR)/.zig-bin;$$env:PATH"' > .env.ps1

.env.sh:
	@echo "export PATH=$(CURDIR)/.zig-bin:\$$PATH" > .env.sh

# COMPILER INSTALLATION TARGETS ================================================================================================================================

ifeq ($(OS),windows)
.zig-bin/zig.exe: | .zig-bin
	@curl.exe -L# -o zig.zip https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).zip ; tar -xf zig.zip -C .zig-bin --strip-components=1 ; rm zig.zip
else
.zig-bin/zig: | .zig-bin
	@curl -L# https://ziglang.org/download/$(ZIG_VERSION)/zig-$(ARCH)-$(OS)-$(ZIG_VERSION).tar.xz | tar -Jx -C .zig-bin --strip-components=1
endif

# DIRECTORY CREATION TARGETS ===================================================================================================================================

.zig-bin:
	@$(if $(filter windows,$(OS)),mkdir .zig-bin -Force | Out-Null,mkdir -p .zig-bin)

example/basis:
	@$(if $(filter windows,$(OS)),mkdir example/basis -Force | Out-Null,mkdir -p example/basis)

# ADDITIONAL TARGETS ===========================================================================================================================================

clean:
	@git clean -dffx

# COMPUTE TARGETS ==============================================================================================================================================

ifeq ($(OS),windows)
molecule: zinq
	@Get-ChildItem example/molecule/*.xyz | ForEach-Object { ./zig-out/bin/zinq.exe hf $$_.FullName --direct --optimize ; Move-Item -Force ($$_.FullName + ".opt") $$_.FullName }
else
molecule: zinq
	@for file in example/molecule/*.xyz; do ./zig-out/bin/zinq hf "$$file" --direct --optimize && mv "$$file.opt" "$$file"; done
endif

# BASIS TARGETS ================================================================================================================================================

BASIS ?= \
    sto-2g \
    sto-3g \
    sto-4g \
    sto-5g \
    sto-6g \
    3-21g \
    6-21g \
    6-31g \
    6-31g* \
    6-31g** \
    6-31+g \
    6-31+g* \
    6-31+g** \
    6-31++g \
    6-31++g* \
    6-31++g** \
    6-311g \
    6-311g* \
    6-311g** \
    6-311+g \
    6-311+g* \
    6-311+g** \
    6-311++g \
    6-311++g* \
    6-311++g** \
    def2-svp \
    def2-svpd \
    def2-tzvp \
    def2-tzvpd \
    def2-tzvpp \
    def2-tzvppd \
    def2-qzvp \
    def2-qzvpd \
    def2-qzvpp \
    def2-qzvppd \
    cc-pvdz \
    cc-pvtz \
    cc-pvqz \
    cc-pv5z \
    aug-cc-pvdz \
    aug-cc-pvtz \
    aug-cc-pvqz \
    aug-cc-pv5z \
    def2-svp-rifit \
    def2-svpd-rifit \
    def2-tzvp-rifit \
    def2-tzvpd-rifit \
    def2-tzvpp-rifit \
    def2-tzvppd-rifit \
    def2-qzvp-rifit \
    def2-qzvpp-rifit \
    def2-qzvppd-rifit \
    cc-pvdz-rifit \
    cc-pvtz-rifit \
    cc-pvqz-rifit \
    cc-pv5z-rifit \
    aug-cc-pvdz-rifit \
    aug-cc-pvtz-rifit \
    aug-cc-pvqz-rifit \
    aug-cc-pv5z-rifit

ifeq ($(OS),windows)
basis: | example/basis
	@'$(BASIS)'.Split() | ForEach-Object { bse get-basis $$_ gaussian94 > "example/basis/$$($$_ -replace '\*','s' -replace '\+','p').g94" }
else
basis: | example/basis
	@set -f; for b in $(BASIS); do bse get-basis "$$b" gaussian94 > "example/basis/$$(echo "$$b" | sed 's/\*/s/g; s/+/p/g').g94"; done
endif

# PROFILE TARGETS ==============================================================================================================================================

ifeq ($(OS),linux)
callgrind: DEBUG := 1
profile:   DEBUG := 1

callgrind: zinq
	@valgrind -q --tool=callgrind --callgrind-out-file=callgrind.data ./zig-out/bin/zinq $(ARGS)

profile: callgrind
	@gprof2dot -e 1 -f callgrind -n 5 $(if $(strip $(NODE)$(ROOT)),-z '$(or $(NODE),$(ROOT))') callgrind.data | dot -T svg -o profile.svg
endif
