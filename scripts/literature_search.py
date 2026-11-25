"""Documented literature search supporting the gap statement.

Queries the OpenAlex works endpoint (open, no key required) with the
query strings recorded here, stores the hit counts and the first page of
results for each query, and writes ``results/summary/literature_search.csv``
and a Markdown record used in the supporting information. Web of Science
and Scopus were not available; the query strings are written so that they
can be rerun on those services.
"""

from __future__ import annotations

import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "summary"
USER_AGENT = "launch-cumulative-screening literature search"

QUERIES = {
    "Q1_cumulative_launch_noise": (
        '"cumulative" AND ("launch noise" OR "rocket noise") AND ("frequency" OR "cadence")'
    ),
    "Q2_point_process_disturbance": (
        '("point process" OR "Poisson process") AND ("disturbance" AND "wildlife") AND ("recovery time" OR "return time")'
    ),
    "Q3_coverage_process_disturbance": (
        '("coverage process" OR "busy period" OR "M/G/infinity") AND ("disturbance" OR "noise events")'
    ),
    "Q4_metric_scaling_cumulative_effects": (
        '("cumulative effects" OR "cumulative impact") AND ("number of events" OR "event frequency") AND ("NEPA" OR "environmental assessment") AND (noise OR disturbance)'
    ),
    "Q5_launch_site_environmental_assessment_cadence": (
        '("launch cadence" OR "launch rate" OR "launches per year") AND ("environmental impact" OR "environmental assessment") AND (noise OR "sonic boom")'
    ),
    "Q6_sobol_cumulative_noise": (
        '("Sobol" OR "global sensitivity") AND ("cumulative noise" OR "day-night average sound level")'
    ),
    "Q7_sonic_boom_cumulative_cdnl": ('("sonic boom" AND ("C-weighted day-night" OR "CDNL") AND cumulative)'),
    "Q8_repeated_disturbance_frequency_shorebirds": (
        '(shorebird OR wader OR plover) AND ("disturbance frequency" OR "disturbance rate") AND (foraging OR "return")'
    ),
}


def search(query: str, per_page: int = 25) -> dict:
    params = {"search": query, "per-page": per_page, "sort": "relevance_score:desc"}
    url = "https://api.openalex.org/works?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.load(resp)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(UTC).strftime("%Y-%m-%d")
    rows = []
    md = [
        f"# Literature search record ({stamp})",
        "",
        "Service: OpenAlex works endpoint (api.openalex.org), relevance sort, 25 results inspected per query.",
        "",
    ]
    for key, q in QUERIES.items():
        try:
            res = search(q)
        except Exception as exc:  # noqa: BLE001
            print(f"{key}: FAILED {exc}")
            rows.append({"query_id": key, "query": q, "date": stamp, "n_hits": "", "status": f"failed: {exc}"})
            continue
        n = res.get("meta", {}).get("count", 0)
        rows.append({"query_id": key, "query": q, "date": stamp, "n_hits": n, "status": "ok"})
        md.append(f"## {key}\n\nQuery: `{q}`\n\nHits: {n}\n")
        for w in res.get("results", [])[:10]:
            title = (w.get("title") or "").replace("\n", " ")
            md.append(f"- {w.get('publication_year')}: {title} (DOI: {w.get('doi')})")
        md.append("")
        print(f"{key}: {n} hits")
        time.sleep(0.5)
    with open(OUT_DIR / "literature_search.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=["query_id", "query", "date", "n_hits", "status"])
        wr.writeheader()
        wr.writerows(rows)
    (OUT_DIR / "literature_search.md").write_text("\n".join(md), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
