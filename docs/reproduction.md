# Reproduction

## Environment

Python 3.12 with dependencies locked in `uv.lock`.

```
uv sync --extra dev          # or: make env
```

A `Dockerfile` builds the same environment (`docker build -t launchcum .`).
LaTeX (`latexmk`) and `pandoc` are needed only for the manuscript, supporting
information and comment documents.

## Inputs

```
python scripts/fetch_data.py          # downloads items marked fetch_on_demand and verifies SHA-256
python -c "from launchcum.data import gis; gis.build(pad_lon=-97.15456, pad_lat=25.99630, half_width_km=30.0)"
python -m launchcum.data.digitize     # digitizes the Appendix B contour maps
```

Processed outputs (`data/processed/`) are committed, so the pipeline can run
without network access.

## Pipeline

```
python -m launchcum.pipeline run --config configs/quick.yaml   # about 2 min on a laptop
python -m launchcum.pipeline run --config configs/paper.yaml   # full configuration
```

The full configuration (10,000 realizations per case, Saltelli design with
8,192 base samples, 16 factors and the seed factor at three cadences, 250 m
grid) takes about 15 minutes on an 8-core laptop. Individual stages can be run
with `--only`, for example `--only verification e1_scaling`.

Outputs are written to `results/<name>/` with a JSON sidecar per file holding
the git commit, configuration hash, input checksums, package versions and
seed. `results/summary/` holds the small summaries committed with the
repository, `paper/generated/values.tex` holds every number quoted in the
manuscript and comment.

## Figures, tables and documents

```
make figures                 # paper/figures/fig*.pdf and paper/generated/tables.tex
make paper si comment attachment
```

`make all` runs everything from inputs to documents.

## Checks

```
make test                    # pytest, including closed-form oracles and transcription reconciliation
make lint typecheck          # ruff and mypy --strict on src/
make refs                    # every DOI resolved through Crossref
```

The continuous-integration workflow (`.github/workflows/ci.yml`) runs these
checks and the quick pipeline on every push.
