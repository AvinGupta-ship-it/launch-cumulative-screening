# launch-cumulative-screening

A stochastic-scheduling framework for screening how cumulative-effect metrics
scale with operating frequency, applied to the proposed increase in
Starship-Super Heavy launch and landing operations at the Boca Chica Launch
Site (FAA docket FAA-2026-12808). The repository contains the `launchcum`
package, the transcribed and digitized FAA per-operation inputs, the
experiments, the manuscript and supporting information, and the scoping
comment with its technical attachment.

The framework treats operations as a marked point process on the calendar,
evaluates four classes of metric (additive, energy-averaged, number-above and
exceedance-day, coverage and interval) plus public-access availability, and
quantifies the error of extrapolating per-operation or annual-average
assessments linearly to higher cadence, with a global sensitivity analysis
over the uncertain inputs. Per-event footprints are the FAA's own modelled
values; the results are statements about metric behaviour, not estimates of
actual impacts.

## Installation

```
uv sync --extra dev
```

Python 3.11 or later; dependencies are locked in `uv.lock`. A `Dockerfile` is
provided.

## Reproduction

```
make quick      # reduced configuration, about two minutes
make all        # full configuration, figures, manuscript and comment
```

See `docs/reproduction.md` for the individual steps, `docs/model.md` for the
model specification and `docs/data.md` for data provenance.

## Layout

```
configs/      experiment configurations (YAML, schema-validated)
data/         MANIFEST.yaml, raw inputs, double-entered transcriptions, processed layers
src/launchcum geometry, scheduler, footprints, metrics, closed forms, sensitivity, emissions, pipeline
scripts/      fetch_data.py, check_refs.py, literature_search.py
figures/      one script per manuscript figure and the table generator
tests/        unit, property, validation and regression tests
results/      generated outputs (small summaries committed under results/summary)
paper/        manuscript, supporting information, bibliography, generated values
comment/      scoping comment and technical attachment
```

## Citation

See `CITATION.cff`. Code is released under the BSD 3-Clause License; text,
figures and derived data under CC BY 4.0 (see `LICENSE`).
