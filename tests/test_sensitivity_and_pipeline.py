import re
from pathlib import Path

import numpy as np
import pytest

from launchcum.io import ScenarioConfig, load_pipeline, load_scenario
from launchcum.sensitivity import Factor, ishigami, ishigami_analytic, saltelli_design, sobol_indices

ROOT = Path(__file__).resolve().parents[1]


def test_sobol_on_ishigami_matches_analytic():
    factors = [Factor(f"x{i}", -np.pi, np.pi, 3) for i in range(3)]
    x = saltelli_design(factors, 2048, seed=1)
    y = ishigami(x)
    si = sobol_indices(factors, y, n_boot=50, seed=1)
    ref = ishigami_analytic()
    assert np.allclose(si["S1"].to_numpy(), ref["S1"], atol=0.06)
    assert np.allclose(si["ST"].to_numpy(), ref["ST"], atol=0.06)


def test_scenario_configs_validate_and_label_hypotheticals():
    for p in sorted((ROOT / "configs").glob("scenario_*.yaml")):
        sc = load_scenario(p)
        assert sc.name == p.stem
        if sc.hypothetical:
            assert "hypothetical" in sc.name
            assert "HYPOTHETICAL" in sc.description.upper()


def test_hypothetical_scenario_requires_label():
    with pytest.raises(ValueError):
        ScenarioConfig.model_validate({"name": "scenario_x", "hypothetical": True, "operations": {"n_launches": 10}})


def test_pipeline_configs_load():
    for name in ("paper.yaml", "quick.yaml"):
        cfg = load_pipeline(ROOT / "configs" / name)
        assert cfg.sobol.n_base >= 8
        assert len(cfg.scenarios) == 5
        assert all(Path(p).exists() for p in cfg.scenarios)


def test_manuscript_macros_are_defined():
    values = ROOT / "paper" / "generated" / "values.tex"
    if not values.exists():
        pytest.skip("values.tex not generated; run the pipeline")
    defined = set(re.findall(r"\\newcommand\{\\(val[A-Za-z]+)\}", values.read_text()))
    used = set()
    for tex in list((ROOT / "paper").glob("*.tex")) + list((ROOT / "comment").glob("*.tex")):
        if tex.name == "values.tex":
            continue
        used |= set(re.findall(r"\\(val[A-Za-z]+)", tex.read_text()))
    missing = used - defined
    assert not missing, f"undefined value macros: {sorted(missing)}"
