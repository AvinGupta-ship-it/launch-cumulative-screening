"""Verify every DOI in a BibTeX file against Crossref metadata.

For each entry with a ``doi`` field the Crossref record is fetched and the
title, first author family name and year are compared with the entry.
Entries without a DOI (Federal Register documents, agency reports) are
checked for the presence of the fields required by the citation rules
(volume/page/date/FR document number, or a URL and date). Any mismatch
fails with a non-zero exit code.
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

USER_AGENT = "launch-cumulative-screening reference checker (mailto omitted)"


def parse_bib(text: str) -> list[dict[str, str]]:
    entries = []
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,]+),(.*?)\n\}", text, flags=re.S):
        kind, key, body = m.group(1), m.group(2).strip(), m.group(3)
        fields = {"_type": kind.lower(), "_key": key}
        for fm in re.finditer(r"(\w+)\s*=\s*(\{(?:[^{}]|\{[^{}]*\})*\}|\"[^\"]*\")", body, flags=re.S):
            val = fm.group(2)[1:-1]
            fields[fm.group(1).lower()] = " ".join(val.split())
        entries.append(fields)
    return entries


def norm(s: str) -> str:
    s = re.sub(r"\\[a-zA-Z]+", "", s)
    s = re.sub(r"[{}~]", "", s)
    s = re.sub(r"[^a-z0-9]+", " ", s.lower())
    return s.strip()


def crossref(doi: str) -> dict:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.load(resp)["message"]


def datacite(doi: str) -> dict:
    """DataCite record in the shape of a Crossref message (title, issued, author)."""
    url = "https://api.datacite.org/dois/" + urllib.parse.quote(doi, safe="")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=60) as resp:
        a = json.load(resp)["data"]["attributes"]
    return {
        "title": [t.get("title", "") for t in a.get("titles", [])],
        "issued": {"date-parts": [[a.get("publicationYear")]]},
        "author": [{"family": c.get("familyName") or c.get("name", "")} for c in a.get("creators", [])],
    }


def check_entry(e: dict[str, str]) -> list[str]:
    problems = []
    doi = e.get("doi")
    if doi:
        try:
            msg = crossref(doi)
        except Exception:  # noqa: BLE001
            try:
                msg = datacite(doi)  # Zenodo and other DataCite registrants are not in Crossref
            except Exception as exc:  # noqa: BLE001
                return [f"{e['_key']}: DOI {doi} did not resolve ({exc})"]
        cr_title = norm(msg.get("title", [""])[0]) if msg.get("title") else ""
        bib_title = norm(e.get("title", ""))
        if cr_title and bib_title and not (bib_title[:60] in cr_title or cr_title[:60] in bib_title):
            problems.append(
                f"{e['_key']}: title mismatch\n   bib: {e.get('title')}\n   crossref: {msg.get('title', [''])[0]}"
            )
        yr = None
        for k in ("published-print", "published-online", "issued", "created"):
            parts = msg.get(k, {}).get("date-parts", [[None]])
            if parts and parts[0] and parts[0][0]:
                yr = parts[0][0]
                break
        if yr and e.get("year") and str(yr) != e["year"]:
            # accept online-first vs print year differences of one year
            if abs(int(yr) - int(e["year"])) > 1:
                problems.append(f"{e['_key']}: year {e['year']} vs crossref {yr}")
        authors = msg.get("author") or []
        if authors and e.get("author"):
            fam = norm(authors[0].get("family", ""))
            first_bib = norm(e["author"].split(" and ")[0].split(",")[0])
            if fam and fam not in norm(e["author"]):
                problems.append(f"{e['_key']}: first author {first_bib!r} not found in crossref authors {fam!r}")
        return problems
    t = e["_type"]
    if t == "misc" and e.get("howpublished", "").startswith("Federal Register") or "FR " in e.get("journal", ""):
        for k in ("volume", "pages", "year", "note"):
            if k not in e:
                problems.append(f"{e['_key']}: Federal Register entry missing {k}")
    else:
        if "url" not in e and "doi" not in e:
            problems.append(f"{e['_key']}: no DOI and no URL")
        if "year" not in e:
            problems.append(f"{e['_key']}: missing year")
    return problems


def main(argv: list[str]) -> int:
    path = Path(argv[1]) if len(argv) > 1 else Path("paper/references.bib")
    entries = parse_bib(path.read_text(encoding="utf-8"))
    failures = []
    for e in entries:
        p = check_entry(e)
        status = "ok " if not p else "BAD"
        print(f"{status} {e['_key']}")
        failures += p
        if e.get("doi"):
            time.sleep(0.2)
    for f in failures:
        print(f, file=sys.stderr)
    print(f"{len(entries)} entries, {len(failures)} problems")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
