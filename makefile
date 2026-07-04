DEBUG ?= 0

.PHONY: all zinq run test

all: zinq

zinq:
	@zig build $(if $(filter 0,$(DEBUG)),--release=fast)

run:
	@zig build $(if $(filter 0,$(DEBUG)),--release=fast) run

test:
	@zig build $(if $(filter 0,$(DEBUG)),--release=fast) test

clean:
	@rm -rf .zig-cache build dist lib zig-out zinq.egg-info zinq/bin zigar zigcc zigcpp zigranlib *.json *.mat *.tar.xz
