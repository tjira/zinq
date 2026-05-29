DEBUG ?= 0

all: zinq

zinq:
	@zig build $(if $(filter 0,$(DEBUG)),--release=fast)

run:
	@zig build $(if $(filter 0,$(DEBUG)),--release=fast) run

clean:
	@git clean -dffx
