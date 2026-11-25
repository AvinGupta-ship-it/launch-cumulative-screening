# Mathematical specification of the model

This document specifies the screening model implemented in `launchcum`. It is a
model of stated assumptions whose per-event inputs are the FAA's own modelled
footprints; it does not predict actual impacts at the site.

## 1. Operations as a marked point process

Operations in one scenario year are events $\{(t_i, m_i)\}$ on the calendar
$t \in [0, T)$, $T = 365$ d (366 in leap years), with marks
$m_i = (\text{type}, \text{period}, \text{origin}, \text{pad})$,
type $\in$ {launch, Super Heavy landing, Starship landing, static fire, scrub},
period $\in$ {day, night} with night defined as 22:00–07:00 local as in the
notice of intent (91 FR 65290), origin $\in$ {local, external}.

Three temporal models generate launch times (`launchcum.schedule`):

1. **Homogeneous Poisson.** The annual count is Poisson with mean $N$ and
   times are uniform on the calendar (or on the business-day axis when
   launches are restricted to weekdays). This is the reference case for the
   closed forms in Section 3.
2. **Constrained renewal (weekday).** Exactly $N$ launches are placed on the
   business-day axis (Monday–Friday excluding U.S. federal holidays, resolved
   with the `holidays` package). Each launch is assigned to one of $P$ pads;
   same-pad intervals are $\delta_{\min} + \mathrm{Gamma}(k, \theta)$ with
   $\theta$ chosen so that the expected total occupies a fraction $1 - s$ of
   the axis ($s$ = `renewal_slack`). Interval sets that do not fit the year are
   rejected and redrawn; the rejection count and a feasibility flag are
   returned. Positions on the business-day axis map to calendar days by
   concatenating business days; the model therefore has day resolution, the
   minimum same-pad turnaround is enforced in whole business days, and the
   time of day is drawn independently (Section 1.2).
3. **Campaign (Neyman–Scott).** Parent points arrive as a Poisson process with
   rate $\lambda_c$; each parent has $\mathrm{Poisson}(\bar n)$ offspring
   displaced by $\mathcal N(0, \sigma_c^2)$ days, with
   $\lambda_c \bar n = N$. Offspring are wrapped into the year and, when the
   scenario is weekday-only, generated on the business-day axis.

Each launch spawns a Super Heavy landing at $t_i + \Delta_{SH}$
($\Delta_{SH}$ = 10 min, the interval the 2025 EA states between launch and
booster landing) with probability $p_{SH}$, and a Starship landing on day
$\lfloor t_i + \Delta_{SS} \rfloor$ with $\Delta_{SS} \sim U(\Delta_{\min},
\Delta_{\max})$ and an independent day/night draw, with probability $p_{SS}$;
a landing drawn earlier in the day than its launch is moved to the same time
on the following day.
Landings of vehicles launched elsewhere arrive as an independent Poisson
process with mean $N_{ext}$ on the full calendar. Static fires precede each
launch (Poisson number with mean `static_fires_per_launch`, daytime, 1–14 d
before launch). Scrubbed attempts precede each launch as a geometric number of
events with success probability $1 - p_s$, each an independent 1–3 d
recycle interval before the next attempt. Scrubs and static fires falling before 1 January are dropped.

The time of day of each event is drawn uniformly within its period. The night
indicator is recomputed from the realized time so that Super Heavy landings
inherit the period of the launch they follow.

### Closure windows

Every launch, Starship landing, static fire and scrub carries an
access-restriction window of lognormal duration $h$ with configured mean and
coefficient of variation; Super Heavy landings are covered by the launch
window. The window starts a fraction `lead_fraction` of $h$ before the event.
Each launch independently triggers an anomaly with probability $p_a$ adding a
fixed `anomaly_hours` to the annual total; anomaly hours are not placed on the
calendar and so do not enter the availability metrics. The launch-window mean (12.9 h) and the 5–25 h Sobol
range bracket the per-launch hours stated in the 2025 Tiered EA for Flights
1–4 (24.5, 10.5, 8.75, 7.75 h); the other closure parameters are
hypothetical and are labelled so in the configuration files.

## 2. Per-event footprints

Footprints are taken from the FAA record (`launchcum.footprints`,
`launchcum.data.digitize`):

- Contour maps in Appendix B of the 2025 Tiered EA (KBR TN 23-05) are
  rasterized at 300 dpi; the legend colours are matched to pixels, the scale
  bar sets the pixel scale, and the pad is the centroid of the innermost
  contour. For each contour the median radius over the overland sector
  (azimuth 240°–360°, where the receptors lie) is recorded, together with the
  minimum and maximum over all azimuths.
- Contour extents and receptor-level bands quoted in the EA text are
  transcribed with page citations, double entered and reconciled.
- For each event type and metric ($L_{A\max}$, SEL) the level at distance $r$
  is the monotone interpolation in $\log r$ through the digitized contour
  radii, held at the innermost value inside the innermost contour, and
  extended beyond the outermost contour by the fitted law
  $L(r) = L_0 - 20\log_{10}(r/r_0) - \alpha(r - r_0)$ (shared $\alpha$ per
  metric) shifted to match the outermost contour. Residuals of the law and of
  the hybrid evaluator against every point are written to
  `footprint_fit_residuals.csv`.
- Sonic boom overpressure follows a weighted power law
  $\log_{10} p = a + b\log_{10} r$ through the transcribed contour extents and
  receptor values, capped at the maximum stated in the record (21 psf for
  Super Heavy, 4 psf for Starship landings). A hypothetical azimuth-dependent
  elongation is available but set to zero in all reported runs.
- Hypothetical vehicle variants enter as an additive level offset $s_L$ and a
  multiplicative boom factor $s_p$.

The C-weighted sound exposure level of a boom follows the FAA Order 1050.1F
Desk Reference (chapter 11): $L_{pk} = 127.6 + 20\log_{10}(p_{psf})$,
$\mathrm{CSEL} = L_{pk} - 26$.

## 3. Metric algebra (`launchcum.metrics`, `launchcum.closed_form`)

For a receptor $x$ and a realization:

- Additive: $M = \sum_i q_{m_i}$ with per-operation quantities from the EA
  tables divided by the stated activity counts.
- Energy-averaged: $L_{dn}(x) = 10\log_{10}\big[\tfrac{1}{T_s}\sum_i
  10^{(L_{AE,i}(x) + 10\,\mathbb 1_{night,i})/10}\big]$, $T_s = 365 \times
  86\,400$ s; CDNL uses CSEL in place of $L_{AE}$, which reduces to the Desk
  Reference form $\mathrm{CSEL} + 10\log_{10}(N_d + 10 N_n) - 49.4$ per day.
  Closed form: $\Delta L_{dn} = 10\log_{10}(N_2/N_1) +
  10\log_{10}[(1 + 9 f_{n,2})/(1 + 9 f_{n,1})]$.
- Number-above $NA_L$, night number-above, boom counts $N_p$ for
  $p \in \{1, 2, 4, 6\}$ psf, exceedance days $D_L$ (distinct calendar days with
  an event above $L$). Poisson closed form $E[D_L] = D(1 - e^{-N_L/D})$ with
  $D$ the number of eligible days.
- Coverage: events exceeding the disturbance threshold (90 dBA $L_{A\max}$ or
  1 psf) start intervals $[t_i, t_i + \tau_i]$, $\tau_i$ exponential with mean
  $\tau$; the disturbed fraction of a season window, the longest undisturbed
  gap and the indicator that a gap of length $\ge w$ exists are computed by
  interval union. Closed form $\phi = 1 - e^{-\lambda E[\tau]}$
  ($M/G/\infty$ busy-period argument).
- Additive: deluge water is 422,000 gal per launch, booster landing and
  booster static fire in the 2025 EA; the model applies it to every static
  fire because static fires are not split by vehicle, an upper bound.
- Access: union length of closure windows (hours), anomaly hours, number of
  windows, fraction of weekend-and-holiday day-period hours (07:00–22:00) and
  of summer (Memorial Day–Labor Day) day-period hours not inside a closure
  (exact overlap of the merged windows with each clock hour), and the
  probability that annual hours exceed caps of 500, 800 and 2,000 h. With the
  2022 PEA measures applied, windows overlapping protected holidays or summer
  weekends are removed.

## 4. Aggregation error

For a metric $Y$ and reference cadence $N_0 = 25$,
$\varepsilon_Y(N) = [Y_{lin}(N) - E[Y(N)]]/E[Y(N)]$ with
$Y_{lin}(N) = (N/N_0)\,E[Y(N_0)]$; for $L_{dn}$ the energy-average
extrapolation error in dB is reported instead. Bootstrap intervals over
realizations accompany each error.

## 5. Sensitivity analysis

Sixteen inputs and the random-number seed of the realization are sampled
uniformly over documented ranges with a Saltelli design ($N = 2^{13}$ base
samples, first- and total-order indices, Jansen estimators as implemented in
SALib, 200 bootstrap resamples). Treating the seed as a factor applies the same
random numbers across the Saltelli matrices, so Monte Carlo noise is reported
as the seed's own index instead of inflating the total-order indices of the
model inputs. The model is the
campaign generator on business days at 25, 100 and 456 launches per year; the
outputs are $L_{dn}$, CDNL, number-above, exceedance days and boom counts at
the Port Isabel edge receptor, disturbed fraction and longest gap at the
piping-plover critical-habitat receptor, closure hours and weekend
availability. The implementation is verified on the Ishigami function.

## 6. Fleet stratospheric emissions

Annual emissions above the tropopause for species $s$ are
$E_s = \sum_{sites} N\, M_{prop}\, f_{strat}\, EI_s$ with uniform ranges on
$M_{prop}$ (6,000–6,750 t from the 2025 EA capacities), $f_{strat}$
(0.5–0.8) and $EI_s$ from the cited literature; no atmospheric modelling is
performed.

## 7. Numerical settings

Realizations use independent `numpy.random.Generator` streams spawned from
one `SeedSequence`, so results do not depend on the number of workers. Grid
fields are computed on a 250 m grid over 60 km × 60 km (EPSG:32614) from
expected day and night counts; a 125 m grid provides the convergence check.
