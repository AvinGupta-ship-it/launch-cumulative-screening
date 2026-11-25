"""Double-entry reconciliation of values transcribed from documents.

Every transcribed table exists as two independently entered files,
``<name>_entry_a.csv`` and ``<name>_entry_b.csv``. :func:`reconcile`
compares them row by row on a key and fails on any mismatch, returning the
reconciled table. The reconciled table is the only one used downstream.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

TRANSCRIBED_DIR = Path(__file__).resolve().parents[3] / "data" / "transcribed"

TABLES: dict[str, list[str]] = {
    # table name -> key columns
    "per_operation_values": ["quantity", "event_type"],
    "footprint_points": ["event_type", "metric", "point_id"],
    "closure_record": ["flight"],
    "emission_rates_table3": ["activity", "pollutant"],
    "ghg_table4": ["activity"],
    "scenario_definitions": ["scenario", "quantity"],
    "verification_targets": ["quantity"],
}


class TranscriptionMismatch(ValueError):
    pass


def _load(name: str, entry: str, base: Path) -> pd.DataFrame:
    path = base / f"{name}_entry_{entry}.csv"
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def reconcile(name: str, base: Path = TRANSCRIBED_DIR) -> pd.DataFrame:
    """Return the reconciled table ``name`` or raise :class:`TranscriptionMismatch`."""
    keys = TABLES[name]
    a = _load(name, "a", base)
    b = _load(name, "b", base)
    value_cols = [c for c in a.columns if c not in keys and c in b.columns]
    if set(a.columns) != set(b.columns):
        raise TranscriptionMismatch(f"{name}: column sets differ: {sorted(set(a.columns) ^ set(b.columns))}")
    a = a.sort_values(keys).reset_index(drop=True)
    b = b.sort_values(keys).reset_index(drop=True)
    if len(a) != len(b) or not a[keys].equals(b[keys]):
        raise TranscriptionMismatch(f"{name}: key rows differ")
    problems = []
    for c in value_cols:
        for i, (va, vb) in enumerate(zip(a[c], b[c], strict=True)):
            if va == vb:
                continue
            try:
                if np.isclose(float(va), float(vb), rtol=0, atol=1e-9):
                    continue
            except ValueError:
                pass
            problems.append(f"{name} row {dict(a.loc[i, keys])} column {c!r}: {va!r} != {vb!r}")
    if problems:
        raise TranscriptionMismatch("\n".join(problems))
    return _typed(a)


def _typed(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        col = out[c].replace("", np.nan)
        try:
            out[c] = pd.to_numeric(col)
        except (ValueError, TypeError):
            out[c] = col
    return out


def reconcile_all(base: Path = TRANSCRIBED_DIR) -> dict[str, pd.DataFrame]:
    return {name: reconcile(name, base) for name in TABLES}


def per_operation_dict(df: pd.DataFrame) -> dict[str, dict[str, float]]:
    """Nested mapping ``quantity -> event_type -> value`` from the reconciled table."""
    out: dict[str, dict[str, float]] = {}
    for _, row in df.iterrows():
        out.setdefault(str(row["quantity"]), {})[str(row["event_type"])] = float(row["value"])
    return out
