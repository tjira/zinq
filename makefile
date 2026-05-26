all: build

bootstrap:
	@julia -e 'using Pkg; Pkg.add(["BenchmarkTools", "PackageCompiler", "Revise"])'

build: instantiate
	@julia --project=. -e 'using PackageCompiler; create_app(".", "dist", executables=["zinq" => "julia_main"], include_lazy_artifacts=true, precompile_execution_file="src/Zinq.jl", force=true)'

clean:
	@git clean -dffx

instantiate:
	@julia --project=. -e 'using Pkg; Pkg.instantiate()'

resolve:
	@julia --project=. -e 'using Pkg; rm("Manifest.toml"); Pkg.resolve()'

repl: bootstrap instantiate
	@julia --project=. -e 'using Revise, Zinq, BenchmarkTools' -i

run: instantiate
	@julia --project=. -e 'using Zinq; julia_main()'

update:
	@julia --project=. -e 'using Pkg; Pkg.update()'
