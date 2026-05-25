bootstrap:
	@julia -e 'using Pkg; Pkg.add(["BenchmarkTools", "PackageCompiler", "Revise"])'

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
