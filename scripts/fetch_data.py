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
import urllib.parse
import urllib.request
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
MANIFEST = ROOT / "data" / "MANIFEST.yaml"
USER_AGENT = "Mozilla/5.0 (compatible; launch-cumulative-screening data fetch)"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, query: str | None, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if query and "resultOffset=<k*1000>" in query:
        _fetch_paginated(url, query, dest)
        return
    full = url if not query else f"{url}?{_encode(query)}"
    req = urllib.request.Request(full, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=600) as resp, open(dest, "wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)


def _encode(query: str) -> str:
    pairs = [kv.split("=", 1) for kv in query.split("&")]
    return urllib.parse.urlencode({k: v for k, v in pairs})


def _fetch_paginated(url: str, query: str, dest: Path) -> None:
    feats = []
    offset = 0
    while True:
        q = query.replace("<k*1000>", str(offset))
        req = urllib.request.Request(f"{url}?{_encode(q)}", headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=600) as resp:
            j = json.load(resp)
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
            # Live services may change; report rather than silently accept.
            print(f"{'CHANGED':9s} {rel}: sha256 {got} differs from manifest {want}")
            failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
