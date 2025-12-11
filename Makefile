ifndef PYTHON_VERSION
PYTHON_VERSION=3.11
endif

sync:
	uv sync

lstm.md: lstm.qmd
	uv run quarto render lstm.qmd --to gfm --output lstm.md

install_python:
	uv python install $(PYTHON_VERSION) && uv sync --python $(PYTHON_VERSION)
