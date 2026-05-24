instantiate:
	julia --project=. -e "using Pkg; Pkg.instantiate()"

resolve:
	julia --project=. -e "using Pkg; rm(\"Manifest.toml\"); Pkg.resolve()"
