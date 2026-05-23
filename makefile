.PHONY: lint

lint: fix format pyright mypy

fix:
	@ruff check --fix .

format:
	@ruff format .

pyright:
	@pyright .

mypy:
	@mypy .
