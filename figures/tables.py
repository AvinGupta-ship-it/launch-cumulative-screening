"""Generate the LaTeX tables used by the manuscript and the comment attachment.

Writes ``paper/generated/tables.tex`` with one macro-free ``tabular`` per
table so the same file can be ``\\input`` by the manuscript, the
supporting information and the technical attachment.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import COMMUNITY_EDGE, COMMUNITY_SPI, HABITAT, RECEPTOR_LETTERS, ROOT, read_csv, read_parquet, summary

OUT = ROOT / "paper" / "generated" / "tables.tex"

EVENT_LABELS = {
    "launch": "Launch",
    "ss_landing": "Starship landing",
    "sh_landing": "Super Heavy landing",
    "static_fire": "Static fire (Super Heavy)",
    "static_fire_ss": "Static fire (Starship)",
}
METRIC_LABELS = {"lamax": r"$L_{A\max}$", "lmax": r"$L_{\max}$ (unweighted)", "sel": "SEL"}


def row(*cells: object) -> str:
    """Join table cells into one LaTeX row."""
    return " & ".join(str(c) for c in cells) + r" \\"


def esc(s: object) -> str:
    return str(s).replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")


def fmtnum(v: float, nd: int = 2) -> str:
    if abs(v) >= 1000:
        return f"{v:,.0f}"
    if abs(v) >= 100:
        return f"{v:.0f}"
    if abs(v) >= 10:
        return f"{v:.1f}"
    return f"{v:.{nd}f}"


def math(s: str) -> str:
    """Set a signed number in math mode so that the minus sign is a true minus."""
    return f"${s}$"


def table_scenarios() -> str:
    rows = [
        (
            "Existing (2025 Tiered EA)",
            "25",
            "25",
            "25",
            r"$\leq 3$ night launches; night landings offshore",
            "25 launches, 3 at night; Poisson",
        ),
        ("Phase 1", "100", "100", "100", "any day or time; 50\\% night assumed", "100 launches; Poisson; 50\\% night"),
        (
            "Phase 2",
            "456",
            "456",
            "456",
            "Monday--Friday only; 50\\% night assumed",
            "456 launches; weekday renewal, 2 pads (hypothetical)",
        ),
        (
            "Phase 2, clustered",
            "456",
            "456",
            "456",
            "as Phase 2",
            "Neyman--Scott campaigns on business days (hypothetical)",
        ),
        (
            "Phase 2, hypothetical upscaled vehicle",
            "456",
            "456",
            "456",
            "as Phase 2",
            "+6 dB levels, 1.5$\\times$ boom (hypothetical)",
        ),
    ]
    out = [
        r"\begin{tabular}{p{2.6cm}rrrp{3.3cm}p{3.6cm}}",
        r"\hline",
        row(
            "Scenario",
            "Launches/yr",
            "SH landings/yr",
            "SS landings/yr",
            "Timing stated in the record",
            "Model configuration",
        ),
        r"\hline",
    ]
    for r in rows:
        out.append(" & ".join(r) + r" \\")
    out += [r"\hline", r"\end{tabular}"]
    return "\n".join(out)


def table_per_operation() -> str:
    po = read_csv("per_operation_inputs.csv")
    fits = read_csv("footprint_fits.csv")
    types = ["launch", "sh_landing", "ss_landing", "static_fire"]
    quantities = [
        ("deluge_gal", "Deluge water (gal per operation)"),
        ("co_tons", "CO (short tons per operation)"),
        ("nox_tons", r"NO$_x$ (short tons per operation)"),
        ("co2e_t", r"CO$_2$e (metric tons per operation)"),
    ]
    out = [
        r"\begin{tabular}{P{3.4cm}rrrr}",
        r"\hline",
        row("Quantity", "Launch", "SH landing", "SS landing", "Static fire"),
        r"\hline",
    ]
    for q, label in quantities:
        cells = []
        for et in types:
            v = po[(po["quantity"] == q) & (po["event_type"] == et)]["value"]
            cells.append(fmtnum(float(v.iloc[0])) if len(v) else "--")
        out.append(row(label, *cells))
    out += [
        r"\hline",
        r"\multicolumn{5}{p{13.2cm}}{Sources. Deluge water: 2025 Tiered EA p.~39 (422,000 gal per launch, "
        r"booster landing or booster static fire; applied to every static fire, an upper bound). CO and NO$_x$: "
        r"2025 Tiered EA Table 3 annual totals divided by 25 operations (11 static fires, the Appendix B scenario). "
        r"CO$_2$e: 2025 Tiered EA Table 4 annual totals divided by the same counts; the launch value includes "
        r"methane venting.} \\",
        r"\hline",
        r"\multicolumn{5}{p{13.2cm}}{Footprint fits. Digitized contour radii are interpolated exactly; the law "
        r"$L(r)=L_0-\alpha\,(r-1\,\mathrm{km})$ continues them beyond the outermost contour.} \\",
        row("Event type and metric", "$L_0$ at 1 km (dB)", r"$\alpha$ (dB/km)", "Points", "rms (dB)"),
        r"\hline",
    ]
    for _, r in fits.iterrows():
        out.append(
            row(
                f"{EVENT_LABELS.get(str(r['event_type']), esc(r['event_type']))}, {METRIC_LABELS[str(r['metric'])]}",
                f"{r['l0_db_at_1km']:.1f}",
                f"{r['alpha_db_per_km']:.2f}",
                int(r["n_points"]),
                f"{r['rms_residual_db']:.1f}",
            )
        )
    out += [
        r"\hline",
        r"\multicolumn{5}{p{13.2cm}}{$\alpha$ is one shared slope per metric, fitted jointly across event types. "
        r"$L_{\max}$ (unweighted) is the unweighted maximum sound pressure level used in the record for structural "
        r"screening. The rms is the root-mean-square difference between the footprint and all transcribed points: "
        r"0.0 dB means every point is a digitized contour radius, which the interpolant reproduces exactly; nonzero "
        r"values arise where text statements at named communities (for example 100--110 dBA) accompany the digitized "
        r"radii and the interpolant passes between them.} \\",
        r"\end{tabular}",
    ]
    return "\n".join(out)


def table_receptors_named() -> str:
    rec = read_csv("receptors_used.csv").set_index("name")
    out = [
        r"\begin{tabular}{clrrP{4.2cm}}",
        r"\hline",
        row("Key", "Receptor", "Distance (km)", "Bearing (deg)", "How located"),
        r"\hline",
    ]
    for letter, name in RECEPTOR_LETTERS.items():
        r = rec.loc[name]
        out.append(
            row(letter, esc(name), f"{float(r['r_km']):.1f}", f"{float(r['azimuth_deg']):.0f}", esc(r["derivation"]))
        )
    out += [
        r"\hline",
        r"\multicolumn{5}{p{12.3cm}}{Bearing is measured clockwise from north at the pad reference point "
        r"(OpenStreetMap apron centroid). Place boundaries are TIGER/Line 2025; the critical-habitat polygon is "
        r"the \mbox{USFWS} designation. Letters key the receptors in Figure~2.} \\",
        r"\end{tabular}",
    ]
    return "\n".join(out)


def table_metric_definitions() -> str:
    rows = [
        (
            "Additive",
            r"$M=\sum_i q_{m_i}$",
            r"$E[M]=\sum_\tau N_\tau \bar q_\tau$",
            "exact; linear",
        ),
        (
            "Energy-averaged",
            r"$\mathrm{DNL}=10\log_{10}\big[\tfrac{1}{T}\sum_i E_i\big]$,\newline "
            r"$E_i=10^{(L_{AE,i}+10\,\mathbb{1}_{n,i})/10}$",
            r"$\Delta=10\log_{10}(N_2/N_1)$ $+\,10\log_{10}\tfrac{1+9f_{n,2}}{1+9f_{n,1}}$",
            "logarithmic in count",
        ),
        (
            "Number-above",
            r"$NA_L=\#\{i: L_{A\max,i}>L\}$",
            r"$E[NA_L]=\sum_\tau N_\tau\,\mathbb{1}[L_\tau>L]$",
            "linear",
        ),
        (
            "Exceedance days",
            r"$D_L=\#\{\text{days with an event}>L\}$",
            r"$E[D_L]=D\,(1-e^{-N_L/D})$",
            "saturates at eligible days $D$",
        ),
        ("Boom count", r"$N_p=\#\{i: p_i>p\}$", "linear", "linear"),
        (
            "Coverage",
            r"$\phi=$ fraction of season in $\bigcup_i[t_i,t_i+\tau_i]$",
            r"$\phi=1-e^{-\lambda E[\tau]}$",
            r"saturates; linear only if $\lambda E[\tau]\ll 1$",
        ),
        (
            "Interval",
            r"longest gap; $P(\text{no gap}\geq w)$",
            r"$\approx\exp\{-\lambda S e^{-\lambda(w+\tau)}\}$",
            "depends on clustering",
        ),
        (
            "Access",
            r"$H=$ union of closure windows; $P(H>H_{cap})$",
            r"$E[H]\approx N\,\bar h$ minus overlap",
            "linear minus overlaps",
        ),
    ]
    out = [
        r"\begin{tabular}{P{1.9cm}P{4.6cm}P{4.8cm}P{2.3cm}}",
        r"\hline",
        r"Class & Definition & Poisson closed form & Scaling with count \\",
        r"\hline",
    ]
    for r in rows:
        out.append(" & ".join(r) + r" \\")
    out += [
        r"\hline",
        r"\multicolumn{4}{p{14.2cm}}{$q_{m_i}$: per-operation quantity of event $i$; $N_\tau$: number of events of "
        r"type $\tau$; $f_n$: night share; $L_{AE}$: sound exposure level; $\mathbb{1}_{n,i}$: night indicator; "
        r"$\lambda$: event rate; $\tau$: recovery time; $S$: season length; $w$: window length; $\bar h$: mean "
        r"closure window.} \\",
        r"\end{tabular}",
    ]
    return "\n".join(out)


def table_verification() -> str:
    vt = read_csv("verification.csv")
    out = [
        r"\begin{tabular}{P{4.3cm}rrP{4.6cm}}",
        r"\hline",
        r"Quantity & 2025 EA value & Model & Basis \\",
        r"\hline",
    ]
    names = {
        "annual NOX (short tons)": "Annual NO$_x$ (short tons)",
        "annual CO (short tons)": "Annual CO (short tons)",
        "annual CO2E (metric tons)": "Annual CO$_2$e (metric tons)",
        "deluge water per operation (gallons)": "Deluge water per operation (gallons)",
    }
    for _, r in vt.iterrows():
        ea_raw = str(r["ea_value"])
        ea = ea_raw.replace("-", "--")
        rep = fmtnum(float(r["reproduced"]))
        try:
            ea_num = float(ea_raw)
            decimals = len(ea_raw.split(".")[1]) if "." in ea_raw else 0
            ea = f"{ea_num:,.{decimals}f}"
            rv = float(r["reproduced"])
            if decimals:
                rep = f"{rv:,.{decimals}f}"
            elif abs(rv - round(rv)) < 1e-9:
                rep = f"{rv:,.0f}"
            else:
                rep = fmtnum(rv)
        except ValueError:
            pass
        out.append(f"{names.get(str(r['quantity']), esc(r['quantity']))} & {ea} & {rep} & {esc(r['basis'])} \\\\")
    out += [r"\hline", r"\end{tabular}"]
    return "\n".join(out)


def _ldn_energy_rel_error(top: float, ref: float) -> tuple[float, float, float]:
    """Relative error of the energy-average extrapolation in the energy domain, with a bootstrap interval."""
    col = f"ldn|{COMMUNITY_EDGE}"
    raw = read_parquet("e1_raw.parquet", columns=["cadence", "external_fraction", col])
    raw = raw[raw["external_fraction"] == 0.0]
    a = raw[raw["cadence"] == ref][col].to_numpy(float)
    b = raw[raw["cadence"] == top][col].to_numpy(float)
    rng = np.random.default_rng(0)

    def rel(ya: float, yb: float) -> float:
        return float(10 ** ((ya + 10 * np.log10(top / ref) - yb) / 10) - 1)

    point = rel(a.mean(), b.mean())
    boots = [rel(rng.choice(a, a.size).mean(), rng.choice(b, b.size).mean()) for _ in range(400)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return point, float(lo), float(hi)


def table_aggregation_errors() -> str:
    err = read_csv("e1_aggregation_errors.csv")
    ci = read_csv("e1_aggregation_error_intervals.csv")
    top = err["cadence"].max()
    e = err[err["cadence"] == top].merge(ci[ci["cadence"] == top], on=["metric", "cadence"], how="left")
    labels = {
        "deluge_gal": "Deluge water (million gal/yr)",
        "nox_tons": "NO$_x$ (short tons/yr)",
        f"ldn|{COMMUNITY_EDGE}": "DNL, Port Isabel edge (dB)",
        f"na100|{COMMUNITY_EDGE}": r"Events $>100$ dBA $L_{A\max}$, Port Isabel edge (per yr)",
        f"days90|{COMMUNITY_EDGE}": r"Days $>90$ dBA, Port Isabel edge (days/yr)",
        f"days90|{COMMUNITY_SPI}": r"Days $>90$ dBA, S. Padre Island edge (days/yr)",
        f"boom2|{COMMUNITY_EDGE}": r"Booms $>2$ psf, Port Isabel edge (per yr)",
        f"disturbed|{HABITAT}": "Disturbed fraction of season, habitat, 1 h recovery",
        "closure_h": "Closure hours, operations (h/yr)",
    }
    out = [
        r"\begin{tabular}{P{4.6cm}rrrr}",
        r"\hline",
        f"Metric at {int(top)}/yr & Simulated & Extrapolated from 25/yr & Error & "
        r"Rel. error, \% (95\% CI) \\",
        r"\hline",
    ]
    ref = float(err["cadence"].min())
    for _, r in e.iterrows():
        if r["metric"] not in labels:
            continue
        scale = 1e-6 if r["metric"] == "deluge_gal" else 1.0
        sim = float(r["simulated"]) * scale
        lin = float(r["linear_extrapolation"]) * scale
        er = float(r["error"]) * scale
        if abs(sim) < 1:
            f_val, f_err = "{:.3f}", "{:+.3f}"
        elif abs(sim) < 100:
            f_val, f_err = "{:.1f}", "{:+.2f}"
        else:
            f_val, f_err = "{:,.0f}", "{:+,.1f}"
        if r["metric"].startswith("ldn"):
            pt, lo, hi = _ldn_energy_rel_error(float(top), ref)
            rel = f"{math(f'{100 * pt:+.1f}')} ({math(f'{100 * lo:+.1f}')}, {math(f'{100 * hi:+.1f}')})"
        else:
            rel = math(f"{100 * r['relative_error']:+.1f}")
            if not pd.isna(r.get("rel_error_lo", float("nan"))):
                lo_s = math(f"{100 * r['rel_error_lo']:+.1f}")
                hi_s = math(f"{100 * r['rel_error_hi']:+.1f}")
                rel += f" ({lo_s}, {hi_s})"
        out.append(
            row(
                labels[r["metric"]],
                f_val.format(sim),
                f_val.format(lin),
                math(f_err.format(er).replace(",", "{,}")),
                rel,
            )
        )
    out += [
        r"\hline",
        r"\multicolumn{5}{p{13.6cm}}{Error is the extrapolated value minus the simulated value; relative error is "
        r"the error divided by the simulated value. Extrapolation is linear in the count from the 25/yr simulation, "
        r"or $+10\log_{10}(456/25)$ for DNL, whose relative error is given for the energy "
        r"$10^{\mathrm{DNL}/10}$. Intervals are bootstrap 95\% intervals over the simulated years. The two "
        r"exceedance-day rows are identical because every launch and landing exceeds 90 dBA at both receptors "
        r"and static fires exceed it at neither, so the same event days are counted; their intervals differ only by "
        r"bootstrap resampling.} \\",
        r"\end{tabular}",
    ]
    return "\n".join(out)


def table_receptors() -> str:
    rec = read_csv("receptors_used.csv")
    out = [r"\begin{tabular}{lrrl}", r"\hline", r"Receptor & $r$ (km) & azimuth (deg) & kind \\", r"\hline"]
    for _, r in rec.iterrows():
        out.append(f"{esc(r['name'])} & {float(r['r_km']):.1f} & {float(r['azimuth_deg']):.0f} & {esc(r['kind'])} \\\\")
    out += [r"\hline", r"\end{tabular}"]
    return "\n".join(out)


def table_factor_ranges() -> str:
    f = read_csv("e6_factor_ranges.csv")
    out = [r"\begin{tabular}{lrrrp{7cm}}", r"\hline", r"Factor & low & high & rank & note \\", r"\hline"]
    for _, r in f.iterrows():
        out.append(
            row(esc(r["name"]), f"{float(r['low']):g}", f"{float(r['high']):g}", int(r["source_rank"]), esc(r["note"]))
        )
    out += [r"\hline", r"\end{tabular}"]
    return "\n".join(out)


def table_temporal_models() -> str:
    d = read_csv("e2_differences.csv")
    top = d["cadence"].max()
    d = d[d["cadence"] == top]
    keep = [
        f"days90|{COMMUNITY_SPI}",
        f"disturbed|{HABITAT}",
        f"longest_gap|{HABITAT}",
        f"has_window|{HABITAT}",
        "weekend_avail",
        "summer_avail",
        "closure_h",
    ]
    labels = {
        f"days90|{COMMUNITY_SPI}": "Days $>90$ dBA, S. Padre Island edge",
        f"disturbed|{HABITAT}": "Disturbed fraction, habitat (1 h recovery)",
        f"longest_gap|{HABITAT}": "Longest gap, habitat (d)",
        f"has_window|{HABITAT}": "P(4 h window), habitat",
        "weekend_avail": "Weekend daylight available",
        "summer_avail": "Summer daylight available",
        "closure_h": "Closure hours (operations)",
    }
    out = [
        r"\begin{tabular}{p{4cm}lrrr}",
        r"\hline",
        f"Metric at {int(top)} launches/yr & Model & Mean & Difference vs Poisson & 95\\% interval \\\\",
        r"\hline",
    ]
    for m in keep:
        for _, r in d[d["metric"] == m].iterrows():
            ci = f"({fmtnum(float(r['ci_lo']), 3)}, {fmtnum(float(r['ci_hi']), 3)})"
            out.append(
                row(labels[m], esc(r["model"]), fmtnum(float(r["mean"]), 3), fmtnum(float(r["diff_vs_poisson"]), 3), ci)
            )
    out += [r"\hline", r"\end{tabular}"]
    return "\n".join(out)


def table_cap_probabilities() -> str:
    t = read_csv("e4_cap_probabilities.csv")
    out = [
        r"\begin{tabular}{p{2.8cm}rrrrrr}",
        r"\hline",
        r"Distribution & Cadence & Mean h & $P(H>500)$ & $P(H>800)$ & $P(H>2000)$ & Weekend avail. \\",
        r"\hline",
    ]
    for _, r in t.iterrows():
        out.append(
            row(
                esc(str(r["distribution"]).replace("_", " ")),
                int(r["cadence"]),
                fmtnum(float(r["closure_anom_h_mean"])),
                f"{float(r['p_exceed_500']):.3f}",
                f"{float(r['p_exceed_800']):.3f}",
                f"{float(r['p_exceed_2000']):.3f}",
                f"{float(r['weekend_avail']):.3f}",
            )
        )
    out += [r"\hline", r"\end{tabular}"]
    return "\n".join(out)


def table_dominant_factors() -> str:
    d = read_csv("e6_dominant_factors.csv")
    out = [r"\begin{tabular}{rllrr}", r"\hline", r"Cadence & Output & Factor & $S_T$ & 95\% half-width \\", r"\hline"]
    for _, r in d.iterrows():
        out.append(
            row(
                int(r["cadence"]),
                esc(r["output"]),
                esc(r["factor"]),
                f"{float(r['ST']):.2f}",
                f"{float(r['ST_conf']):.2f}",
            )
        )
    out += [r"\hline", r"\end{tabular}"]
    return "\n".join(out)


def table_emissions() -> str:
    e = read_csv("e8_emissions.csv")
    out = [
        r"\begin{tabular}{llrrrr}",
        r"\hline",
        r"Site & Species & Launches/yr & p05 (t) & median (t) & p95 (t) \\",
        r"\hline",
    ]
    for _, r in e.iterrows():
        out.append(
            row(
                esc(r["site"]),
                esc(r["species"]),
                int(r["launches_per_year"]),
                fmtnum(float(r["p05_t"])),
                fmtnum(float(r["median_t"])),
                fmtnum(float(r["p95_t"])),
            )
        )
    out += [r"\hline", r"\end{tabular}"]
    return "\n".join(out)


def main() -> None:
    s = summary()
    parts = [
        "% Generated by figures/tables.py from saved results; do not edit by hand.",
        r"\providecommand{\tabularnewline}{\\}",
        r"\newcolumntype{P}[1]{>{\raggedright\arraybackslash}p{#1}}",
        r"\newcommand{\tableScenarios}{" + table_scenarios() + "}",
        r"\newcommand{\tablePerOperation}{" + table_per_operation() + "}",
        r"\newcommand{\tableReceptorsNamed}{" + table_receptors_named() + "}",
        r"\newcommand{\tableMetricDefinitions}{" + table_metric_definitions() + "}",
        r"\newcommand{\tableVerification}{" + table_verification() + "}",
        r"\newcommand{\tableAggregationErrors}{" + table_aggregation_errors() + "}",
        r"\newcommand{\tableReceptors}{" + table_receptors() + "}",
        r"\newcommand{\tableFactorRanges}{" + table_factor_ranges() + "}",
        r"\newcommand{\tableTemporalModels}{" + table_temporal_models() + "}",
        r"\newcommand{\tableCapProbabilities}{" + table_cap_probabilities() + "}",
        r"\newcommand{\tableDominantFactors}{" + table_dominant_factors() + "}",
        r"\newcommand{\tableEmissions}{" + table_emissions() + "}",
        f"% pipeline {s.get('pipeline_name')} seed {s.get('seed')}",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(parts) + "\n", encoding="utf-8")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
