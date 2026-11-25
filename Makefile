PY ?= .venv/bin/python
UV ?= uv

.PHONY: all quick env data digitize test lint typecheck refs figures paper si comment attachment clean

all: env data digitize
	$(PY) -m launchcum.pipeline run --config configs/paper.yaml
	$(MAKE) figures
	$(MAKE) paper si comment attachment

quick: env data digitize
	$(PY) -m launchcum.pipeline run --config configs/quick.yaml

env:
	$(UV) sync --extra dev

data:
	$(PY) scripts/fetch_data.py
	$(PY) -c "from launchcum.data import gis; gis.build(pad_lon=-97.15456, pad_lat=25.99630, half_width_km=30.0)"

digitize:
	$(PY) -m launchcum.data.digitize > /dev/null

test:
	$(PY) -m pytest -q

lint:
	$(PY) -m ruff check src tests scripts figures
	$(PY) -m ruff format --check src tests scripts figures

typecheck:
	$(PY) -m mypy --strict src

refs:
	$(PY) scripts/check_refs.py paper/references.bib

figures:
	for f in figures/fig*.py; do $(PY) $$f; done

paper:
	cd paper && latexmk -pdf -interaction=nonstopmode manuscript.tex
	cd paper && pandoc manuscript.tex --citeproc --bibliography=references.bib -o manuscript.docx

si:
	cd paper && latexmk -pdf -interaction=nonstopmode si.tex

comment:
	cd comment && latexmk -pdf -interaction=nonstopmode comment.tex
	cd comment && pandoc comment.tex -o comment.docx

attachment:
	cd comment && latexmk -pdf -interaction=nonstopmode attachment_technical.tex

clean:
	rm -rf results/quick results/paper
	cd paper && latexmk -C manuscript.tex si.tex || true
	cd comment && latexmk -C comment.tex attachment_technical.tex || true
