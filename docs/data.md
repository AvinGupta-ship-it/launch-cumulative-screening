# Data sources, variables and limitations

All inputs are listed with retrieval date, URL, query and SHA-256 in
`data/MANIFEST.yaml`. Values read from documents are stored in
`data/transcribed/` as two independently entered files per table and are
reconciled by `launchcum.data.transcribed.reconcile` (a test fails on any
mismatch). Every transcribed row cites the document page.

## Source hierarchy

1. Statute, regulation, Federal Register documents and docket documents
   (notice of intent, Federal Register notices of availability, the 2022 PEA,
   the 2025 Tiered EA and its appendices, the 2026 Tiered EA, the FAA scoping
   presentation, the FAA Desk Reference).
2. Agency datasets with documented methodology (U.S. Census TIGER/Line and
   TIGERweb 2020 blocks, USFWS critical-habitat and refuge-boundary services,
   USGS PAD-US, TPWD state-park boundaries served by TxDOT).
3. Peer-reviewed literature (emission indices, disturbance ecology, point
   processes, sensitivity analysis).
4. Government technical reports (FICON 1992).
5. Grey literature and crowd-sourced data (the launch-pad reference point from
   OpenStreetMap). Any rank-5 quantity is flagged and included in the
   sensitivity analysis or shown to be immaterial.

## Regulatory record (rank 1)

| Document | Used for | Limitation |
|---|---|---|
| 91 FR 65289 (NOI, FR Doc 2026-20963) and 91 FR 72611 (update, FR Doc 2026-22898) | Cadence scenarios (25/100/456 launches; 912 landings regardless of origin), weekday restriction, 50% night assumption, day/night definition, resource areas, alternatives, deadline | No per-operation quantities |
| FAA scoping presentation (November 2026) | 2,000 h/yr access-restriction assumption; landing locations | Presentation, not a decision document |
| 2025 Final Tiered EA and FONSI/ROD | Operations table (22 day/3 night), emission tables 3 and 4, deluge water, closure caps and the Flight 1–4 closure record, noise contour statements | Modelled, operator-supplied inputs |
| 2025 Tiered EA Appendix B (KBR TN 23-05; sonic boom analysis) | Contour maps digitized for $L_{A\max}$ and SEL of launches, landings and static fires; DNL map for the 25-operation scenario; CDNL basis (44 booms) | Contours are for the vehicle configuration modeled in 2024 |
| 2026 Tiered EA (trajectories, Starship RTLS) | Confirms caps unchanged (500 + 300 h); no new quantitative footprints used | — |
| 91 FR 5020 and 90 FR 59094 | Fleet cadence for the emissions inventory (44 and 76 launches/yr) | — |
| FAA Desk Reference ch. 11 | psf to CSEL/CDNL conversion | — |

Inconsistencies found in the record and how they are handled:

- The 2025 EA text states the 100 dBA SEL contour extends about 8 miles west
  and the 90 dBA contour about 15 miles; the digitized Appendix B maps give
  overland radii of about 14.8 and 27 miles. The maps are the modeled
  output and are used; the text extents are retained at reduced weight and
  their residuals are reported.
- Structural-damage contours are labeled "dBA" in places; they are treated
  as unweighted $L_{\max}$ with the other structural contours and are not used
  for A-weighted metrics.

## Spatial data (rank 2)

- Census 2020 blocks with population (TIGERweb) for population counts and the
  Boca Chica Village receptor; TIGER/Line 2025 places for community receptors
  (internal points and nearest polygon edges).
- USFWS critical habitat: piping plover final designation (74 FR 23476) and
  rufa red knot proposed designation (88 FR 22530; status "Proposed" in the
  service at retrieval).
- USFWS refuge boundaries (Lower Rio Grande Valley and Laguna Atascosa NWR).
- Boca Chica State Park boundary (TPWD provisional boundary via TxDOT) and
  PAD-US units for Isla Blanca Park.
- TIGER/Line 2025 area hydrography for Cameron County (Laguna Madre, Gulf,
  South Bay, Bahia Grande), used only for the shoreline outlines on the maps.
- Launch pad reference point: OpenStreetMap apron centroid; a 300 m
  uncertainty changes the nearest community distance by under 4%.

## Literature values (rank 3)

- Recovery times: Collop et al. (2016) report 27–75 s mean (up to 163 s) to
  resume feeding after pedestrian disturbance for ten wader species; van der
  Kolk et al. (2020) report GPS-tracked oystercatchers returning "often
  within a few hours" after heavy aircraft disturbance; Stewart (1993)
  reports pinnipeds returning within several hours after Titan IV launches.
  The sweep covers 1 min to 2 days; no species-specific value is asserted for
  the site.
- Disturbance-frequency effects: Lafferty (2001), Tarr et al. (2010), van der
  Kolk et al. (2024), Weston and Elgar (2007).
- Emission indices: Barker et al. (2024) methane row; Ross and Sheaffer
  (2014) and Ryan et al. (2022) for brackets and the stratospheric fraction;
  Tsigaridis et al. (2024) for the methane black-carbon assumption; Maloney et
  al. (2022) for comparison scenarios.

## Hypothetical parameters

The following enter only as labeled hypothetical ranges: same-pad turnaround,
number of pads, campaign clustering parameters, scrub probability, static
fires per launch, Starship landing delay, closure window parameters other
than the launch mean, anomaly probability, vehicle-variant scaling factors,
boom directivity, the required undisturbed window $w$, and the share of
Starship landings at the site.

## Not available

A complete official record of closure hours could not be compiled: Cameron
County's notice archive is no longer served, archived notices give scheduled
windows rather than realized hours, and no agency publishes a tally. Closure
durations are therefore parameterized from the per-flight hours stated in the
2025 EA and the gap is stated in the manuscript and the comment.
