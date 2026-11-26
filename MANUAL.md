---
title: "Frequency Scaling of Cumulative Effects for High-Cadence Launch and Landing Operations"
subtitle: "Project Technical Design Manual"
date: "Regulatory status verified November 23, 2026"
---

# 1. Project Specification

This manual is the complete technical specification for one computational research project. Its finished form is a reproducible research repository, a scientific preprint, and a formal regulatory comment. Section 1 specifies the project (A: identity and significance; B: methodology; C: data and references; D: codebase; E: experiments; F: validation and limitations; G: preprint; H: regulatory comment; I: final artifacts). Section 2 sets the scientific standards the work must meet. Section 3 defines the deliverable architecture.

**Conventions.** Every statement of regulatory status was verified from the Federal Register, Regulations.gov, agency web pages or docket documents on November 23, 2026, and must be re-verified on the day any comment is submitted. Sources marked **[verify]** were identified but not confirmed (DOI not resolved, document not read, or URL not tested); they must be confirmed before use, and may not be cited until they are. Facts drawn from secondary sources are labeled as such. Numbers quoted from Federal Register documents were read from text extractions and must be checked against the official PDF before they are quoted in a manuscript or comment. All scenario parameters not taken from a cited source are hypothetical and are labeled as such in code, figures and text.


## 1.A Project identity and research significance

**Manuscript title.** *How impact metrics scale with operating frequency: a stochastic-scheduling framework for cumulative-effects screening of high-cadence launch and landing operations*

**Agency and docket.** Federal Aviation Administration, Office of Commercial Space Transportation. Docket FAA-2026-12808.

**Verified proceeding.** Notice of Intent to Prepare a Programmatic Environmental Impact Statement (PEIS), Open Public Scoping Period, and Hold Public Scoping Meetings on the Proposed SpaceX Starship-Super Heavy Launch and Landing Operations at the Boca Chica Launch Site at Starbase, Texas and Expansion of Starship Landing and Expenditure Locations. 91 FR 65289–65291 (October 14, 2026), FR Doc. 2026-20963, Regulations.gov document FAA-2026-12808-0001; updated by 91 FR 72611–72612 (November 10, 2026), FR Doc. 2026-22898, document FAA-2026-12808-0057, which extended scoping and converted the meetings to virtual format.

**Status (verified November 23, 2026).** The scoping comment period is **open until December 7, 2026** (postmark or electronic submission). No Draft PEIS has been released and no date for one has been published. The scoping meetings (November 17–19, 2026) have occurred. A Draft PEIS comment period will follow at a date not yet announced; that is the second, more quantitative opportunity in the same proceeding.

**Federal action under review.** Modification of SpaceX's vehicle operator license (14 CFR part 400) for Starship–Super Heavy at Boca Chica, subsequent modifications and renewals, and FAA-issued temporary airspace closures, reviewed under NEPA and FAA Order 1050.1G. The proposed action, from the NOI and the FAA scoping presentation:

| Scenario | Launches/yr | Super Heavy landings/yr | Starship landings/yr | Timing |
|---|---|---|---|---|
| Existing (2025 Tiered EA) | 25 | 25 | 25 | ≤ 3 night launches; night landings offshore only |
| Phase 1 | 100 | 100 | 100 | any day or time |
| Phase 2 ("within several years") | 456 | 456 | 456 | any time of day, Monday–Friday only |

Additional elements: up to 912 landings/yr at Boca Chica "regardless of launch origin" (vehicles launched from other sites may land at Starbase); noise modeling assumes 50% of operations at night (22:00–07:00); public access restrictions of up to **2,000 hours/yr** analyzed, against the 2025 limit of 500 hours/yr for operations plus 300 hours/yr for anomalies; new Pacific and Indian Ocean landing and expenditure areas. Alternatives are the Proposed Action and No Action only. Resource areas listed for scoping include air quality, biological resources, coastal resources, DOT Act §4(f) and LWCF §6(f), hazardous materials, cultural resources, land use, noise and compatible land use, socioeconomics and children's health, visual effects, and water resources. Cooperating agencies: NASA, NPS, USCG, USFWS.

**Exact regulatory question.** What environmental variables, metrics, alternatives and analyses should the PEIS use to characterize the effects of increasing operations from 25 to 100 to 456 launches per year (and up to 912 landings), given that the existing NEPA record for the site is built largely on per-operation footprints and annual totals?

**Scientific problem.** Environmental assessments of launch sites have been built around single-event footprints (maximum sound level, sonic boom overpressure contours, per-launch emissions, per-operation water use and closure duration) and annual sums or energy averages of those footprints. Whether that construction remains adequate when frequency increases by more than an order of magnitude depends on how each effect aggregates over events. Some effects are additive (emissions mass, water volume); some are logarithmic in event count (energy-averaged noise metrics); some saturate (fraction of time a receptor is disturbed); some depend on the time between events (recovery of wildlife behavior, availability of a public beach on a given day), and these respond to clustering, not only to the count. A screening framework that makes this aggregation structure explicit is generic to high-cadence operations and does not presently exist in a form tied to launch-site NEPA metrics **[novelty to be confirmed by literature review; see 1.G]**.

**Central research question.** For each class of cumulative-effect metric, how does the metric scale with annual operating frequency and temporal clustering, how large is the error from extrapolating per-operation assessments linearly or by energy averaging, and which uncertain inputs control the cumulative result at each cadence?

**Original contribution.**

1. A taxonomy of cumulative-effect metrics by aggregation operator (additive, energy-averaged, coverage/saturating, interval/threshold), with closed-form scaling laws under Poisson operations and corrections for clustered (weekday, campaign) operations.
2. A stochastic operations scheduler that produces calendar-resolved realizations of launches, landings, static fires, scrubs and closures consistent with stated constraints (weekday-only, day/night split, pad turnaround), from which any metric can be computed with uncertainty.
3. A quantitative demonstration, using the FAA's own published per-operation values as inputs, of where per-operation and annual-average metrics diverge from interval-based metrics as cadence increases, and a global sensitivity analysis ranking the inputs that control each metric.

**Independent value.** The aggregation results apply to any repeated-disturbance source (airports with growing operations, offshore pile driving, military ranges, wind-turbine construction). The scheduler and metric library are reusable for other launch sites, including the Florida Starship sites.

**Relationship to the agency's decision.** At scoping, the FAA decides what the PEIS will analyze. The project's findings bear directly on which metrics, sensitivity analyses and alternatives are needed to characterize the proposed cadence; they do not and cannot establish the magnitude of real impacts, which requires the Draft PEIS's modeling.

**Chain of relevance.** Research question (how metrics scale with cadence) → stochastic scheduling plus metric algebra using FAA per-operation inputs → finding on which metrics change behavior with cadence and which inputs dominate → scoping issue (metrics, analyses and alternatives for the PEIS) → recommendation on PEIS metrics, disclosures and alternatives.

## 1.B Mathematical and computational methodology

**B.1 Operations as a marked point process.** Let operations on the calendar \(t\in[0,T]\), \(T=365\) d, be events \(\{(t_i, m_i)\}\) with marks \(m_i=(\text{type}, \text{period}, \text{origin})\), type ∈ {launch, Super Heavy landing, Starship landing, static fire, scrub}, period ∈ {day, night}. Three stochastic models of operations are implemented and compared:

1. *Homogeneous Poisson* with rate \(\lambda=N/T\): the reference case for closed forms.
2. *Constrained renewal (weekday) model*: launches occur only on Monday–Friday business days \(\mathcal{D}\); inter-launch intervals on the same pad follow a shifted gamma distribution with minimum turnaround \(\delta_{\min}\) and shape \(k\); number of pads \(P\in\{1,2\}\) is a parameter. Each launch spawns a Super Heavy landing at \(t_i+\Delta_{\text{SH}}\) (minutes) and a Starship landing at \(t_i+\Delta_{\text{SS}}\) (hours to days, or none for expendable missions). Landings of vehicles launched elsewhere arrive as an independent process with rate \(\lambda_{\text{ext}}\) up to the 912/yr total.
3. *Campaign (clustered) model*: a Neyman–Scott cluster process with cluster rate \(\lambda_c\), mean cluster size \(\bar n\), and intra-cluster spread \(\sigma_c\), matched to the same annual count, to represent campaign-style operations such as propellant-transfer sequences.

Scrubs and recycles are Bernoulli trials per attempt with probability \(p_s\), each consuming a closure window without a launch. Night fraction \(f_n=0.5\) per the NOI is the baseline; \(f_n\in[0,0.5]\) is varied.

**B.2 Per-event footprints.** Footprints are taken from the FAA record, not from a new acoustic or atmospheric model:

- *Launch noise* \(L_{A\max}(r)\) and sound exposure level \(L_{AE}(r)\) as functions of distance from the pad, digitized from contour maps and tables in the 2022 PEA and the 2025 Tiered EA, interpolated with a fitted form \(L(r)=L_0-20\log_{10}(r/r_0)-\alpha (r-r_0)\), with residuals reported. Static-fire levels likewise.
- *Sonic boom* peak overpressure for Super Heavy landings from the 2025 Tiered EA (up to 21 psf near the site; approximately 6 psf at 10 mi, 4 psf at 15 mi, 2 psf at 27 mi, 1 psf at 28 mi, as transcribed in 1.C) fitted by a piecewise log-linear interpolant in \((\log r,\log p)\); Starship landing booms 1–4 psf.
- *Per-operation quantities*: deluge water (≈422,000 gal per launch, landing or static fire in the 2025 EA), annual criteria and greenhouse-gas emissions at 25 launches (NOx 98 tpy, CO 45 tpy, 97,342 t CO2e), converted to per-operation values with the operation mix stated.
- *Closure duration*: the distribution of public-access restriction hours per operation, scrub and anomaly, estimated from public closure notices if a complete record can be compiled (1.C); otherwise a documented parameter range bracketing the 2025 cap divided by planned operations.

The model is explicitly a screening model: footprints are FAA-modeled values for the current vehicle configuration, and their extension to future configurations is a stated assumption tested by scaling factors.

**B.3 Metric algebra.** For a receptor at position \(x\) and a realization of operations, four metric classes are computed.

*(i) Additive (mass/volume)*: \(M=\sum_i q_{m_i}\). Exact scaling \(E[M]=\sum_{\text{types}} N_\tau \bar q_\tau\); only uncertainty in \(\bar q_\tau\) matters.

*(ii) Energy-averaged noise*: community day–night level
\[
L_{dn}(x)=10\log_{10}\!\left[\frac{1}{T_s}\sum_i 10^{(L_{AE,i}(x)+10\,\mathbb{1}_{\text{night},i})/10}\right],\qquad T_s=86{,}400\ \text{s}\times 365 ,
\]
and its C-weighted analogue (CDNL) for launch noise as used in the FAA record. Under fixed per-event levels, \(L_{dn}\) increases by \(10\log_{10}(N_2/N_1)\) plus a night-weighting term \(10\log_{10}[(1+9f_{n,2})/(1+9f_{n,1})]\). The 25→456 increase with \(f_n\) moving from ≈0.12 to 0.5 implies about 12.6 dB from count and roughly 3–4 dB from the night share; the code reproduces these closed forms as tests and maps the resulting contour expansion.

*(iii) Number-above and exceedance-day metrics*: \(NA_L(x)=\#\{i: L_{A\max,i}(x)>L\}\); \(D_L(x)=\#\{\text{days with at least one event above }L\}\); boom counts \(N_{p}(x)=\#\{i:p_i(x)>p\}\) for \(p\in\{1,2,4,6\}\) psf. Under Poisson operations, \(E[D_L]=365\,(1-e^{-\lambda_L})\) with \(\lambda_L\) the daily rate of exceeding events; under weekday operations the bound is the number of weekdays, so the per-day exposure intensity rather than the day count grows at high cadence. Both are reported.

*(iv) Coverage and interval metrics for disturbance-sensitive receptors*: each event disturbs a receptor for a recovery time \(\tau\) (random, distribution \(G\)). The fraction of time disturbed is the coverage of the union of intervals \([t_i,t_i+\tau_i]\); for a Poisson process (an \(M/G/\infty\) busy-period argument) the stationary probability of being in a disturbed state is
\[
\phi = 1-\exp(-\lambda\,E[\tau]),
\]
which is linear in \(\lambda\) only when \(\lambda E[\tau]\ll 1\). The undisturbed-interval distribution, the probability that no undisturbed window of length \(w\) occurs within a season (relevant to nesting or foraging bouts), and the longest undisturbed gap are computed by simulation and compared with the Poisson closed forms. For clustered processes the same mean rate gives larger undisturbed gaps but more intense disturbed periods; the framework quantifies which effect dominates for given \(\tau\) and \(w\). Recovery times are not known for the species and stressors at the site; \(\tau\) and \(w\) are swept over ranges reported in the disturbance literature (1.C) and the results are presented as response surfaces, not species-specific predictions.

*(v) Public access availability*: from the closure process, compute annual closure hours \(H\), the probability \(P(H>H_{\text{cap}})\) for caps of 500, 800 and 2,000 h, the fraction of weekend and holiday daylight hours available (the Phase 2 weekday restriction matters here), and the availability during the summer peak season. The ratio \(H_{\text{cap}}/N\) gives the mean closure budget per launch implied by the cap; the model reports the probability the cap is consistent with a cadence given the per-operation closure distribution.

**B.4 Aggregation-error estimator.** For a metric \(Y\), define the linear-extrapolation error relative to a reference cadence \(N_0\):
\[
\varepsilon_Y(N)=\frac{Y_{\text{lin}}(N)-E[Y(N)]}{E[Y(N)]},\qquad Y_{\text{lin}}(N)=\frac{N}{N_0}E[Y(N_0)] ,
\]
and the analogous error for an energy-average extrapolation. These errors are the central quantitative object of the paper.

**B.5 Global sensitivity analysis.** Inputs \(\boldsymbol\theta\) (night fraction, scrub probability, turnaround, number of pads, clustering parameters, footprint scaling factors \(s_L, s_p\) for future vehicle variants, closure duration parameters, \(\tau\), \(w\), fraction of landings from other sites) are sampled from documented ranges. First-order and total Sobol indices are estimated for each metric at each cadence (Saltelli sampling, \(N=2^{13}\) base samples; bootstrap intervals). The output is a cadence-by-metric map of dominant uncertainties.

**B.6 Stratospheric emissions (bounded sub-analysis).** For the combined Starship cadence across permitted and proposed sites (Starbase Phase 2, 456/yr; KSC LC-39A, 44/yr per the January 2026 ROD; SLC-37, 76/yr per the November 2025 Department of the Air Force ROD), compute annual emissions of H2O, CO2, NOx and black carbon above the tropopause from methane–oxygen combustion using emission indices and altitude partitioning from the literature, with ranges, and compare the order of magnitude to the scenarios in Maloney et al. (2022) and Ryan et al. (2022). This is restricted to an emissions inventory and comparison; no new atmospheric chemistry modeling is performed.

**Numerical methods.** Event simulation is vectorized in NumPy; each configuration uses \(R=10^4\) annual realizations (or more until the Monte Carlo standard error criterion of Section 2.2 is met). Spatial metrics are computed on a 250 m grid over a 60 km × 60 km domain in UTM zone 14N (EPSG:32614), with receptor points for named communities (Port Isabel, South Padre Island, Brownsville, Boca Chica Village), Boca Chica State Park, and critical-habitat polygons.

## 1.C Data and scientific references

**Docket and NEPA record (primary).**

| Source | Identifier | Content used | Entry into model | Limitations |
|---|---|---|---|---|
| NOI, FAA-2026-12808 | 91 FR 65289; FR Doc 2026-20963; https://www.federalregister.gov/documents/2026/10/14/2026-20963/ | Cadence scenarios, night share, landings from other sites, resource areas, alternatives | Scenario definitions | No per-operation quantities |
| NOI update | 91 FR 72611; FR Doc 2026-22898 | Deadline December 7, 2026; virtual meetings | Submission eligibility | — |
| FAA scoping presentation | https://www.faa.gov/space/stakeholder_engagement/spacex_starship/Starbase-PEIS-Virtual-Public-Scoping-Meeting-PPT.pdf | 2,000 h/yr access restrictions; trajectories; azimuths | Closure cap scenarios | Presentation, not a decision document |
| FAA project page | https://www.faa.gov/space/stakeholder_engagement/spacex_starship | Document index | — | — |
| 2022 Final PEA and Mitigated FONSI | 87 FR 36199 (June 15, 2022), FR Doc 2022-12888; https://www.faa.gov/sites/faa.gov/files/2022-06/PEA_for_SpaceX_Starship_Super_Heavy_at_Boca_Chica_FINAL.pdf | Baseline footprints, noise contours, anomaly hours | Footprint digitization | Earlier vehicle configuration |
| 2025 Increased Cadence Tiered EA and FONSI/ROD | 90 FR 20203 (May 12, 2025), FR Doc 2025-08232; EA https://www.faa.gov/media/94346; FONSI/ROD https://www.faa.gov/media/94336 | 25 launches/50 landings; 500+300 h closures; LAmax 100–140 dBA to ≈8 mi; 60 dB CDNL ≈5 mi; SH booms to 21 psf; 422,000 gal per operation; NOx 98 tpy, CO 45 tpy, 97,342 t CO2e | Per-operation inputs and validation targets | Modeled, operator-supplied inputs |
| Tiered EA, Additional Trajectories and Starship RTLS, FONSI/ROD | FR Doc 2026-03291 (February 19, 2026); EA and FONSI PDFs on FAA project page | Landing trajectories, booms from Starship return | Footprint inputs | — |
| KSC LC-39A Final EIS and ROD | FR Doc 2026-02108 (February 3, 2026); https://www.faa.gov/space/stakeholder_engagement/spacex_starship_ksc | 44 launches/yr | B.6 fleet cadence | — |
| SLC-37 ROD (Department of the Air Force) | FR Doc 2025-23314 (December 18, 2025) | 76 launches, 152 landings/yr | B.6 fleet cadence | Revised ROD pending FAA SEIS (docket FAA-2026-10792, scoping closed October 14, 2026) |

Per-operation values listed above were read from the 2025 FONSI/ROD and must be re-transcribed from the source tables (page and table cited) into `data/transcribed/` with double entry before use.

**Spatial and environmental data.**

- USFWS Environmental Conservation Online System critical-habitat data (piping plover wintering units on the lower Texas coast; rufa red knot) — https://ecos.fws.gov **[verify exact layer URLs and current designations]**.
- Lower Rio Grande Valley National Wildlife Refuge and Boca Chica State Park boundaries — USFWS cadastral data and Texas Parks and Wildlife GIS **[verify]**.
- U.S. Census TIGER/Line places and 2020 block population for receptor population counts — https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html.
- Cameron County public notices of State Highway 4 and Boca Chica Beach closures **[verify archive availability and completeness]**; if a complete multi-year record cannot be compiled from official notices, closure durations are parameterized and the gap is stated.
- Texas Open Beaches Act, Tex. Nat. Res. Code ch. 61, and the Cameron County–Texas General Land Office beach-closure memorandum of agreement referenced in the 2025 EA **[verify document]**.

**Scientific literature.**

- Maloney, C. M., Portmann, R. W., Ross, M. N., & Rosenlof, K. H. (2022). The climate and ozone impacts of black carbon emissions from global rocket launches. *Journal of Geophysical Research: Atmospheres*, 127, e2021JD036373. https://doi.org/10.1029/2021JD036373
- Ryan, R. G., Marais, E. A., Balhatchet, C. J., & Eastham, S. D. (2022). Impact of rocket launch and space debris air pollutant emissions on stratospheric ozone and global climate. *Earth's Future*, 10, e2021EF002612. https://doi.org/10.1029/2021EF002612 (metadata consistent with ADS and PMC records; resolve DOI before citing)
- Gee, K. L., et al. Measurements of Starship Super Heavy launch noise and booster-catch events, *JASA Express Letters* (2024–2025) **[verify title and DOI]**; Gee et al. (2023) on Space Launch System sound power, *JASA Express Letters* 3(11) **[verify DOI]**.
- Ross, M. N., & Sheaffer, P. M. (2014). Radiative forcing caused by rocket engine emissions. *Earth's Future*, 2, 177–196 **[verify DOI]**.
- Federal Interagency Committee on Noise (1992), *Federal Agency Review of Selected Airport Noise Analysis Issues* **[verify URL]** — basis of DNL significance conventions.
- Literature on wildlife disturbance and recovery (shorebird flushing and return times; disturbance-frequency effects) **[to be compiled; each value used for τ and w must come from a located, read source]**.
- Daley–Vere-Jones, *An Introduction to the Theory of Point Processes* (Springer) and standard \(M/G/\infty\) results for coverage processes — textbook references for B.3(iv).

## 1.D Research codebase specification

Repository `launch-cumulative-screening`, package `launchcum`.

| Module | Function |
|---|---|
| `launchcum.schedule` | Poisson, constrained-renewal and Neyman–Scott generators; business-day calendar (U.S. federal holidays via `holidays`); scrubs; external-origin landings; returns a typed `pandas.DataFrame` of events |
| `launchcum.footprints` | Loads transcribed FAA footprint tables and digitized contours; fits and evaluates \(L_{A\max}(r)\), \(L_{AE}(r)\), \(p(r)\); azimuth-dependent boom footprints for landing trajectories |
| `launchcum.metrics` | Additive totals; \(L_{dn}\)/CDNL; number-above; exceedance days; boom counts; coverage, gap and longest-gap statistics; access-availability metrics |
| `launchcum.closed_form` | Poisson closed forms used as test oracles and in the manuscript |
| `launchcum.spatial` | Grid and receptor geometry (GeoPandas, Shapely, PyProj); habitat-polygon aggregation |
| `launchcum.sensitivity` | SALib wrappers, Sobol and Morris designs, bootstrap |
| `launchcum.emissions` | Fleet stratospheric emission inventory (B.6) |
| `launchcum.io` | Configuration schemas (pydantic), result writers (Parquet, GeoPackage, NetCDF via xarray) |

Libraries: NumPy, SciPy, pandas, xarray, GeoPandas, Shapely, PyProj, rasterio, SALib, matplotlib, cartopy (maps), pydantic, pint, pytest, hypothesis.

Inputs: YAML configurations (`configs/scenario_existing.yaml`, `phase1.yaml`, `phase2.yaml`, `phase2_clustered.yaml`, `phase2_hypothetical_upscaled_vehicle.yaml`), transcribed tables (CSV with page/table citation columns), vector layers (GeoPackage). Outputs: per-realization metric tables (Parquet), gridded metric fields (NetCDF), receptor summaries (CSV), figures (PDF/PNG).

Tests: closed-form \(L_{dn}\) increments; Poisson coverage \(1-e^{-\lambda E\tau}\) recovered within Monte Carlo error; weekday generator never places launches on weekends or holidays in Phase 2; annual counts match configuration exactly (or in distribution, per model); footprint interpolants reproduce transcribed points within stated tolerance; Sobol indices on the Ishigami function match analytical values.

## 1.E Computational experiments

| # | Varied | Output | Purpose | Evaluation |
|---|---|---|---|---|
| E1 | Cadence N ∈ {25, 50, 100, 200, 456} and external landings | All metric classes at receptors | Establish scaling of each metric class | Compare with closed forms; report \(\varepsilon_Y(N)\) with 95% intervals |
| E2 | Temporal model (Poisson, weekday renewal, clustered) at fixed N | Coverage, longest gap, exceedance days, weekend availability | Isolate the effect of clustering at constant count | Differences with bootstrap intervals; test whether clustering changes the ranking of scenarios |
| E3 | Night fraction 0–0.5 | \(L_{dn}\), night number-above | Quantify night-share contribution | Closed-form check; contour area change |
| E4 | Closure duration distribution and caps 500/800/2,000 h | \(P(H>H_{\text{cap}})\), implied per-operation budget, seasonal availability | Test internal consistency of cadence and closure cap | Probability statements with Monte Carlo error |
| E5 | Recovery time τ and required window w | Disturbed fraction, \(P(\text{no window} \ge w)\) | Map conditions under which per-event and cumulative assessments diverge | Response surfaces; identify regimes \(\lambda E[\tau]\gtrsim 0.1\) |
| E6 | All inputs (Sobol) | Each metric at each cadence | Rank uncertainties | Total-order indices with bootstrap intervals; factors with \(S_T>0.05\) reported |
| E7 | Footprint scaling (hypothetical future vehicle variants) | Contour areas, boom counts | Test robustness to vehicle changes | Labeled hypothetical |
| E8 | Fleet cadence across sites | Stratospheric emission totals | Context for global cumulative effects | Order-of-magnitude comparison with literature scenarios |

**Competing hypotheses.** (H1) Energy-averaged and annual-total metrics track the interval and coverage metrics closely enough across the proposed range that additional metrics add little; (H2) they diverge substantially at Phase 2 cadence, so that the choice of metric changes which receptors and scenarios appear most affected; (H3) divergence depends mainly on unknown recovery parameters, in which case the defensible recommendation is a disclosure-and-sensitivity request rather than a specific metric. The analysis determines which holds for each receptor class.

## 1.F Validation and limitations

1. *Numerical correctness*: closed-form oracles (B.3), conservation of event counts, grid-convergence of contour areas (250 m vs 125 m), Monte Carlo convergence.
2. *Model verification*: reproduce the 2025 EA's own reported quantities at 25 launches (60 dB CDNL extent ≈5 mi; annual NOx 98 tpy; deluge totals) from the transcribed per-operation inputs, within transcription and digitization tolerance. A failure indicates a misread input or an aggregation convention mismatch and must be resolved before proceeding.
3. *Empirical validity*: limited. Where public closure notices allow, compare simulated per-operation closure distributions with the observed record. Published far-field noise measurements (Gee et al.) may be compared with the EA-modeled \(L_{A\max}(r)\) as a plausibility check only. No validation of ecological response is possible; ecological outputs are conditional response surfaces.
4. *Regulatory applicability*: findings are framed as recommendations for PEIS scope and metrics. They do not estimate actual impacts of the proposed action.

Principal limitations: dependence on FAA-modeled footprints for the current vehicle; unknown future operational patterns within the stated caps; unknown species-specific recovery parameters; no atmospheric transport modeling.

## 1.G Scientific preprint

**Title.** *How impact metrics scale with operating frequency: a stochastic-scheduling framework for cumulative-effects screening of high-cadence launch and landing operations.*

**Target venues.** *Environmental Impact Assessment Review*; *Journal of the Acoustical Society of America* (if the noise metrics dominate); preprint on EarthArXiv.

**Structure.**

- *Abstract* (≈250 words): problem (per-event assessment under order-of-magnitude cadence growth), method (marked point processes, metric algebra, Sobol analysis), data basis (FAA-published per-operation values for one site), principal quantitative results, limitations.
- *1 Introduction*: cumulative-effects assessment under NEPA; growth of launch cadence; existing launch environmental literature (acoustics, atmospheric emissions); disturbance ecology of repeated stressors; gap statement after a documented search (Web of Science/Scopus query strings recorded in SI).
- *2 Methods*: 2.1 operations models; 2.2 footprints and their provenance; 2.3 metric classes and closed forms; 2.4 aggregation error; 2.5 sensitivity analysis; 2.6 fleet emission inventory.
- *3 Results*: 3.1 verification against the 2025 EA; 3.2 scaling of each metric class with cadence; 3.3 effect of clustering and weekday restriction; 3.4 closure-cap consistency; 3.5 disturbance-coverage regimes; 3.6 dominant uncertainties.
- *4 Discussion*: implications for metric selection in programmatic reviews; generality beyond launch sites; what measurements would most reduce uncertainty.
- *5 Limitations*; *6 Conclusions*.
- *Figures*: (1) site map with receptors, habitat and footprints; (2) schematic of the four aggregation classes; (3) metric vs cadence curves with linear-extrapolation overlays; (4) \(L_{dn}\) contour maps at 25/100/456 with night share; (5) coverage and longest-gap response surfaces in \((\lambda,\tau)\); (6) closure-hours distributions with caps; (7) Sobol total-order heatmap (metric × input × cadence).
- *Tables*: (1) scenarios; (2) per-operation inputs with sources; (3) metric definitions and scaling laws; (4) verification against 2025 EA values; (5) aggregation errors at Phase 2.
- *Supporting information*: derivations of closed forms; digitization method and residuals; convergence studies; full parameter ranges with sources; literature-search protocol.

## 1.H Formal regulatory comment

**Proceeding and route.** Scoping comment on the NOI, docket FAA-2026-12808, via Regulations.gov (https://www.regulations.gov/docket/FAA-2026-12808) or by mail to Ms. Amy Hanson, FAA Environmental Protection Specialist, SpaceX Boca Chica, c/o Kimley-Horn, 767 Eustis Street, Suite 100, St. Paul, MN 55114. Deadline **December 7, 2026**. Comments are posted publicly without change. If the analysis is not complete by that date, no scoping comment is filed on its basis; the work is directed instead at the Draft PEIS comment period in the same docket, whose opening must then be confirmed from the Federal Register.

**Issues within scope.** The NOI requests comment on alternatives, potential impacts and relevant information and analyses for the PEIS. Directly in scope: (a) metrics used for noise and compatible land use, including night operations at 50% share; (b) public access restrictions up to 2,000 h/yr and effects on §4(f)/§6(f) resources and coastal access; (c) biological resources in the adjacent refuge and critical habitat, including cumulative disturbance frequency; (d) water use for deluge at higher cadence; (e) landings of vehicles launched from other sites; (f) the range of alternatives (only Proposed Action and No Action are identified); (g) cumulative effects with other Starship sites.

**Findings that could support recommendations (to be determined by results).**

- Whether energy-averaged noise metrics alone adequately describe the change from 25 to 456 launches, or whether number-above, exceedance-day and boom-count metrics at named receptors change the picture (E1, E3).
- Whether temporal clustering under the weekday restriction materially changes undisturbed-interval statistics for disturbance-sensitive receptors (E2, E5), supporting a request that the PEIS analyze temporal distribution, not just annual counts.
- Whether the 2,000 h cap is internally consistent with the proposed cadence under plausible per-operation closure durations, and what fraction of weekend and peak-season hours would remain available (E4).
- Which uncertain inputs dominate each metric (E6), supporting targeted sensitivity analyses in the PEIS.
- Whether intermediate-cadence alternatives (for example Phase 1 as a separate alternative, or caps on night operations) would be distinguishable on the metrics that matter (E1, E3), supporting a request for additional alternatives.

**Comment structure.** (1) Commenter and docket; (2) one-page summary of recommendations; (3) description of the cadence increase as stated in the NOI; (4) analysis: aggregation classes, scaling results, clustering, closure-cap consistency, sensitivity; (5) recommendations, each tied to a resource area and analysis in the PEIS; (6) limitations, stating that results use the FAA's existing per-operation values and do not estimate actual impacts; (7) references; attachment: technical note with Figures 3–7 and Tables 2–5.

## 1.I Final artifacts

- Repository `launch-cumulative-screening` with tagged release.
- Transcribed FAA per-operation dataset with page-level citations; digitized contours with residuals.
- Simulation outputs for all experiments; gridded metric fields.
- Seven manuscript figures (vector) and five tables generated by scripts.
- Manuscript (LaTeX, PDF, DOCX) and supporting information.
- Scoping comment (DOCX, PDF) and technical attachment, prepared for the December 7, 2026 deadline if results are validated by then; otherwise a Draft-PEIS-stage comment prepared when that period opens.

# 2. Scientific Standards

The standards in this section apply to every part of this project. Where Section 1 is more specific, Section 1 governs.

## 2.1 Computational reproducibility

1. **Single-command reproduction.** Every repository exposes one entry point (`make all`, or `python -m launchcum.pipeline run --config configs/paper.yaml`) that regenerates every number, table and figure in the manuscript and the comment attachments from raw inputs. A second target (`make quick`) runs a reduced configuration (fewer Monte Carlo replicates, coarser grids) for continuous integration and must finish on a laptop.
2. **Pinned environments.** Python ≥ 3.11. Dependencies are declared in `pyproject.toml` and locked (`uv.lock` or `conda-lock.yml`). A `Dockerfile` or `environment.yml` reproduces the environment. Solvers external to Python (HiGHS, AERMOD, AERMET) are pinned by version and checksum, and the build fetches them by script.
3. **Seeded randomness.** All stochastic code takes an explicit `numpy.random.Generator`. Seeds are stored in the configuration file and written into every output's metadata. Parallel Monte Carlo uses `SeedSequence.spawn` so results do not depend on worker count.
4. **Configuration as data.** Every experiment is a YAML file under `configs/` validated by a schema (pydantic). No parameter that appears in the manuscript may be hard-coded in a module.
5. **Output provenance.** Each output file (Parquet, NetCDF, CSV, PNG/PDF) carries or is accompanied by a JSON sidecar recording the git commit, configuration hash, input-data checksums, package versions and seed.
6. **Deterministic figures.** Figures are produced only by scripts in `figures/` reading saved results, never by notebooks with hidden state. Fonts and styles are fixed in a shared `matplotlibrc`.

## 2.2 Numerical validation

Each project separates four levels of evidence, and each manuscript reports them separately:

| Level | Question | Typical evidence |
|---|---|---|
| Numerical correctness | Does the code compute what the equations say? | Unit tests, analytical solutions, conservation and identity checks, convergence under refinement, solver optimality certificates |
| Model verification | Do the equations represent the intended conceptual model? | Comparison with independent implementations or published benchmark cases; limiting-case behavior; dimensional analysis |
| Empirical validity | Does the model reproduce observations of the real system? | Out-of-sample comparison with measured data where any exist; otherwise explicitly stated as not established |
| Regulatory applicability | Does the result bear on the decision the agency is actually making? | Mapping from result to a specific provision, question or criterion in the docket |

Minimum requirements for all projects:

- Monte Carlo estimators report standard errors and the replicate count is chosen so that the Monte Carlo standard error of every reported headline quantity is below 2% of its value (or a stated absolute tolerance), demonstrated by a convergence plot in the supporting information.
- Discretized solvers report a grid- or time-step-refinement study with the observed order of convergence.
- Optimization models report solver status, optimality gap (MILP gap ≤ 0.1% unless stated) and verify primal feasibility of the returned solution with an independent check.
- Units are carried with `pint` in input parsing and checked by dimensional tests at module boundaries; core numerical kernels may operate on plain arrays in SI units.
- Property-based tests (`hypothesis`) cover invariants (non-negativity, monotonicity in a parameter where theory requires it, conservation).

## 2.3 Source reliability and hierarchy

Sources are ranked and the rank is recorded in the bibliography database:

1. Statute, regulation (eCFR), and Federal Register documents; docket documents (technical support documents, supporting statements, memoranda to the docket).
2. Agency datasets with documented methodology (EPA AQS, EIA-930, PHMSA incident data, EMTS, AP-42, USGS).
3. Peer-reviewed literature and recognized technical standards (NFPA, IEC, API, ASTM).
4. National-laboratory and government technical reports (NREL, LBNL, Sandia, GAO, OIG, CRS).
5. Grey literature from trade associations, advocacy organizations and consultancies, used only to locate primary sources or to represent stakeholder positions, never as the sole basis of a quantitative parameter.

A parameter taken from rank-5 material must be flagged in the parameter table and included in the sensitivity analysis.

## 2.4 Data provenance

- Raw downloaded files are stored unmodified under `data/raw/` (or fetched on demand by `scripts/fetch_data.py` when licensing or size forbids redistribution), with a `data/MANIFEST.yaml` recording source organization, URL, retrieval date, SHA-256 checksum, license or terms of use, and the exact query parameters used.
- All transformations from raw to analysis-ready data are code in `src/launchcum/data/` and are idempotent. No manual spreadsheet edits.
- Synthetic data are generated only by code, stored under `data/synthetic/`, and every synthetic file name and metadata record carries the word `synthetic`. Figures based on synthetic data say so in the caption.
- Values digitized or transcribed from documents (tables in a Federal Register notice, a technical support document, a figure in a paper) are stored in `data/transcribed/` with a citation to page and table, and are double-entered and reconciled by a test that compares the two entries.

## 2.5 Citation accuracy

- Every reference is held in a single BibTeX file per repository. Every DOI is resolved through doi.org and the resolved metadata (title, authors, year, venue) is checked against the entry by a script (`scripts/check_refs.py`) that queries Crossref; mismatches fail the build.
- Federal Register documents are cited by volume, page, date and FR document number; regulations by CFR section and the eCFR "current as of" date; docket documents by their Regulations.gov document ID or FERC eLibrary accession number.
- A reference is cited only for a claim it actually supports. Statements in this manual marked **[verify]** identify sources that must be located and checked before use; they may not be cited until verified.
- No preprint DOI, journal acceptance or review status is stated unless it exists.

## 2.6 Statistical rigor

- Uncertainty is propagated, not appended: reported quantities are accompanied by intervals derived from the uncertainty model (parametric distributions, bootstrap, or posterior), with the interval type stated.
- Global sensitivity analysis uses variance-based Sobol indices (Saltelli estimator with Jansen formulas, implemented with `SALib`) when inputs are independent and the model is cheap enough for N(k+2) evaluations with N ≥ 2¹²; otherwise Morris elementary effects for screening followed by Sobol indices on the retained factors. Bootstrap confidence intervals on indices are reported.
- Hypothesis tests are pre-specified in the configuration, with multiplicity controlled (Holm or Benjamini–Hochberg as stated) when more than one comparison supports a single claim.
- Rare-event probabilities below 10⁻³ are estimated with variance reduction (importance sampling or subset simulation) or analytically, never by naive Monte Carlo with too few events. The number of observed events behind each estimate is reported.
- Distributions fitted to data report goodness of fit (Anderson–Darling or a likelihood-based comparison) and the effect of the distributional choice on headline results is tested.

## 2.7 Professional manuscript preparation

- Manuscripts are written in LaTeX (`revtex4-2` or `elsarticle` class, per the intended venue) with a single notation table per paper, SI units, and equations numbered. Every number in the text is generated by the pipeline into `paper/generated/values.tex` as LaTeX macros, so prose and results cannot drift apart. A test fails if any macro is undefined or stale.
- Figures are vector PDF with embedded fonts, colorblind-safe palettes (Okabe–Ito or viridis family), labeled axes with units, and captions that state the data basis (observed, modeled, synthetic, hypothetical scenario).
- Editable output: the LaTeX source plus a DOCX produced with `pandoc` from the same source. PDF output: compiled with `latexmk`.
- Author contributions, funding and competing-interest statements follow venue conventions. Acknowledgments name only people and organizations that contributed.
- Target venues are suggested in each project specification; the manuscript must also stand alone as a preprint (arXiv physics.soc-ph, stat.AP, eess.SY, or EarthArXiv / ESS Open Archive, as appropriate).

## 2.8 Regulatory relevance

Each comment satisfies these tests before it is finalized:

1. **Proceeding test.** The comment is addressed to the exact document and docket identified in the project specification and was confirmed open on the date of submission from the Federal Register or the agency docket. If the period has closed, the comment is not filed as if timely; the specification states the alternative (later stage of the same proceeding, or no filing).
2. **Scope test.** Each recommendation names the provision, question, criterion or decision in the docket it addresses and the agency authority under which it could be acted on in that proceeding.
3. **Evidence test.** Each recommendation is traced to a specific figure or table in the attachment and to the validation level (Section 2.2) that supports it. Recommendations resting only on hypothetical scenarios are phrased as requests for analysis, disclosure or clarification, not as findings about real facilities.
4. **Restraint test.** No claim of harm, violation, or benefit at a real facility is made without facility-specific evidence. Adverse or null results that bear on the question are reported.

Comments follow a common structure: (i) identification of commenter and docket; (ii) summary of recommendations (one page); (iii) the specific issue in the proposal; (iv) analysis and findings, with figures; (v) recommendations, each tied to a provision; (vi) limitations; (vii) references; (viii) technical attachment (a condensed version of the preprint and the repository archive identifier, if public). Comments are submitted as PDF with an editable DOCX retained.

## 2.9 Scientific integrity

- Conclusions are determined by validated results. Each project specification lists competing hypotheses and the observations that would distinguish them; the manuscript reports which outcome occurred, including null or unfavorable results.
- No empirical data, measurements, validation results, citations, or publication metadata are fabricated. A model that cannot be validated against observations is described as a verified model of stated assumptions, not as a prediction.
- Hypothetical scenarios are labeled as such everywhere they appear: configuration names (`scenario_hypothetical_*`), figure captions, table headers, and text.
- Any facility, site, or company is analyzed only from public documents, and inferences about it are stated with the evidence they rest on.

## 2.10 Research repository organization

The repository `launch-cumulative-screening` uses this skeleton:

```
launch-cumulative-screening/
├── README.md               # purpose, installation, reproduction, citation
├── LICENSE                 # BSD-3-Clause for code; CC BY 4.0 for text, figures, derived data
├── CITATION.cff            # citation metadata; version 1.0.0
├── .zenodo.json            # Zenodo deposit metadata (version taken from the release tag)
├── pyproject.toml
├── uv.lock                 # or conda-lock.yml
├── Dockerfile
├── Makefile
├── configs/                # experiment configurations (YAML), schema-validated
├── data/
│   ├── MANIFEST.yaml       # provenance of every input
│   ├── raw/                # unmodified downloads (or fetched by script)
│   ├── transcribed/        # values transcribed from documents, double-entered
│   ├── synthetic/          # generated synthetic inputs
│   └── processed/          # analysis-ready derived data
├── src/launchcum/          # importable library
├── scripts/                # fetch_data.py, check_refs.py, run_experiments.py
├── figures/                # one script per manuscript figure
├── notebooks/              # exploratory analysis with scientific content only
├── tests/                  # unit, property, validation and regression tests
├── results/                # generated outputs (git-ignored except small summaries)
├── docs/
│   ├── model.md            # mathematical specification of the model
│   ├── data.md             # data sources, variables, limitations
│   └── reproduction.md     # step-by-step reproduction and expected runtimes
├── paper/                  # LaTeX manuscript, SI, bibliography, generated values
└── comment/                # regulatory comment source, attachments, final PDF
```

Files that do not serve a scientific, computational, documentation, publication or regulatory function are not created. In particular, repositories contain no development logs, progress records, session notes, change narratives beyond conventional `git` history and a `CHANGELOG.md` limited to released versions, task lists, or files describing how the repository was produced.

Continuous integration (GitHub Actions) runs linting (`ruff`), type checking (`mypy --strict` on `src/`), the test suite, `make quick`, and the reference checker on every push.

# 3. Deliverable Architecture

## 3.1 Repository

The project is implemented as one self-contained repository, `launch-cumulative-screening`, with importable Python package `launchcum`, organized as in Section 2.10 and with the modules listed in Section 1.D. It has no code dependency on any other repository and can be archived, cited and reproduced alone.

## 3.2 Publication package

Each project produces, inside its repository:

- `paper/manuscript.tex`, `paper/manuscript.pdf`, `paper/manuscript.docx` — the preprint.
- `paper/si.tex`, `paper/si.pdf` — supporting information (extended methods, verification results, convergence studies, full parameter tables, additional figures).
- `paper/figures/` — final vector figures as used in the manuscript.
- `comment/comment.docx`, `comment/comment.pdf` — the formal regulatory comment.
- `comment/attachment_technical.pdf` — the technical attachment (condensed methods, results and figures, written for agency technical staff).
- A tagged release (`v1.0.0`) archived to Zenodo when made public, producing a DOI that is then (and only then) cited in the comment and manuscript.
