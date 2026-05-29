DEBUG ?= 0

all: zinq

zinq:
	zig build $(if $(filter 0,$(DEBUG)),--release=fast)

clean:
	git clean -dffx
