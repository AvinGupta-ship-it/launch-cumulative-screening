"""Fetch raw inputs listed in data/MANIFEST.yaml and verify checksums.

Usage::

    python scripts/fetch_data.py            # fetch items marked fetch_on_demand
    python scripts/fetch_data.py --all      # re-fetch everything
    python scripts/fetch_data.py --verify   # only verify checksums of present files
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
MANIFEST = ROOT / "data" / "MANIFEST.yaml"
# Plain client name: the FAA site refuses agent strings with parenthesised comments.
USER_AGENT = "launch-cumulative-screening/1.0"
RETRIES = 4
RETRY_WAIT_S = 65.0


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _get(full: str) -> bytes:
    """GET with retries; ArcGIS error documents (such as rate limits) count as failures."""
    last: Exception | None = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(full, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=600) as resp:
                data = resp.read()
            if data[:1] == b"{" and b'"error"' in data[:200]:
                raise RuntimeError(json.loads(data).get("error", {}).get("message", "service error"))
            return data
        except Exception as exc:  # noqa: BLE001
            last = exc
            if attempt < RETRIES - 1:
                print(f"{'retry':9s} {exc}; waiting {RETRY_WAIT_S:.0f} s")
                time.sleep(RETRY_WAIT_S)
    raise RuntimeError(str(last))


def fetch(url: str, query: str | None, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if query and "resultOffset=<k*1000>" in query:
        _fetch_paginated(url, query, dest)
        return
    full = url if not query else f"{url}?{_encode(query)}"
    dest.write_bytes(_get(full))


def _encode(query: str) -> str:
    pairs = [kv.split("=", 1) for kv in query.split("&")]
    return urllib.parse.urlencode({k: v for k, v in pairs})


def _fetch_paginated(url: str, query: str, dest: Path) -> None:
    feats = []
    offset = 0
    while True:
        q = query.replace("<k*1000>", str(offset))
        j = json.loads(_get(f"{url}?{_encode(q)}"))
        f = j.get("features", [])
        feats += f
        if not j.get("properties", {}).get("exceededTransferLimit") and len(f) < 1000:
            break
        offset += len(f)
    dest.write_text(json.dumps({"type": "FeatureCollection", "features": feats}), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args(argv)
    manifest = yaml.safe_load(MANIFEST.read_text(encoding="utf-8"))
    failures = 0
    changed = 0
    for rel, item in manifest["items"].items():
        dest = RAW / rel
        want = item["sha256"]
        if args.verify:
            if dest.exists():
                got = sha256(dest)
                status = "ok" if got == want else "MISMATCH"
                failures += status != "ok"
                print(f"{status:9s} {rel}")
            else:
                print(f"{'missing':9s} {rel}")
            continue
        if dest.exists() and not args.all and sha256(dest) == want:
            print(f"{'present':9s} {rel}")
            continue
        if not (args.all or item.get("fetch_on_demand")) and dest.exists():
            print(f"{'present':9s} {rel}")
            continue
        print(f"{'fetching':9s} {rel}")
        try:
            fetch(item["url"], item.get("query"), dest)
        except Exception as exc:  # noqa: BLE001
            print(f"{'FAILED':9s} {rel}: {exc}")
            failures += 1
            continue
        got = sha256(dest)
        if got != want:
            # Live services re-serialize or revise their data between releases:
            # report the new checksum so that the manifest can be updated after
            # review, without treating a successful download as a failure.
            print(f"{'CHANGED':9s} {rel}: sha256 {got} differs from manifest {want}")
            changed += 1
    if changed:
        print(f"{changed} file(s) differ from data/MANIFEST.yaml; review and update the manifest if accepted")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
