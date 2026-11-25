"""Generate the LaTeX tables used by the manuscript and the comment attachment.

Writes ``paper/generated/tables.tex`` with one macro-free ``tabular`` per
table so the same file can be ``\\input`` by the manuscript, the
supporting information and the technical attachment.
"""

from __future__ import annotations

import pandas as pd

from common import COMMUNITY_EDGE, COMMUNITY_SPI, HABITAT, ROOT, read_csv, summary

OUT = ROOT / "paper" / "generated" / "tables.tex"


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
        r"\begin{tabular}{lrrrll}",
        r"\hline",
        r"Scenario & Launches/yr & SH landings/yr & SS landings/yr & Timing stated in the record & Model configuration \\",
        r"\hline",
    ]
    for r in rows:
        out.append(" & ".join(r) + r" \\")
    out += [r"\hline", r"\end{tabular}"]
    return "\n".join(out)


def table_per_operation() -> str:
    po = read_csv("per_operation_inputs.csv")
    fits = read_csv("footprint_fits.csv")
    out = [r"\begin{tabular}{llrl}", r"\hline", r"Quantity & Event type & Value & Source \\", r"\hline"]
    src = {
        "deluge_gal": "2025 Tiered EA p.~39 (422,000 gal per operation)",
        "co_tons": "2025 Tiered EA Table 3 / 25 operations",
        "nox_tons": "2025 Tiered EA Table 3 / 25 operations",
        "co2e_t": "2025 Tiered EA Table 4 / 25 operations (launch includes methane venting)",
    }
    for _, r in po.iterrows():
        out.append(
            f"{esc(r['quantity'])} & {esc(r['event_type'])} & {fmtnum(float(r['value']))} & {src.get(r['quantity'], '')} \\\\"
        )
    out += [
        r"\hline",
        r"\multicolumn{4}{l}{Footprint fits (digitized contours, monotone interpolation; law parameters for extrapolation):} \\",
    ]
    for _, r in fits.iterrows():
        out.append(
            f"{esc(r['event_type'])} {esc(r['metric'])} & $L_0$ = {r['l0_db_at_1km']:.1f} dB at 1 km & $\\alpha$ = {r['alpha_db_per_km']:.2f} dB/km & {int(r['n_points'])} points, rms {r['rms_residual_db']:.1f} dB \\\\"
        )
    out += [r"\hline", r"\end{tabular}"]
    return "\n".join(out)


def table_metric_definitions() -> str:
    rows = [
        ("Additive", r"$M=\sum_i q_{m_i}$", r"$E[M]=\sum_\tau N_\tau \bar q_\tau$", "exact; linear"),
        (
            "Energy-averaged",
            r"$L_{dn}=10\log_{10}\!\big[\tfrac{1}{T}\sum_i 10^{(L_{AE,i}+10\,\mathbb{1}_{n,i})/10}\big]$",
            r"$+10\log_{10}(N_2/N_1)+10\log_{10}\tfrac{1+9f_{n,2}}{1+9f_{n,1}}$",
            "logarithmic in count",
        ),
        (
            "Number-above",
            r"$NA_L=\#\{i: L_{A\max,i}>L\}$",
            r"$E[NA_L]=\sum_\tau N_\tau \mathbb{1}[L_\tau>L]$",
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
        r"\begin{tabular}{llll}",
        r"\hline",
        r"Class & Definition & Poisson closed form & Scaling with count \\",
        r"\hline",
    ]
    for r in rows:
        out.append(" & ".join(r) + r" \\")
    out += [r"\hline", r"\end{tabular}"]
    return "\n".join(out)


def table_verification() -> str:
    vt = read_csv("verification.csv")
    out = [r"\begin{tabular}{lrrl}", r"\hline", r"Quantity & 2025 EA value & Reproduced & Basis \\", r"\hline"]
    for _, r in vt.iterrows():
        out.append(
            f"{esc(r['quantity'])} & {esc(r['ea_value'])} & {fmtnum(float(r['reproduced']))} & {esc(r['basis'])} \\\\"
        )
    out += [r"\hline", r"\end{tabular}"]
    return "\n".join(out)


def table_aggregation_errors() -> str:
    err = read_csv("e1_aggregation_errors.csv")
    ci = read_csv("e1_aggregation_error_intervals.csv")
    top = err["cadence"].max()
    e = err[err["cadence"] == top].merge(ci[ci["cadence"] == top], on=["metric", "cadence"], how="left")
    labels = {
        "deluge_gal": "Deluge water (additive)",
        "nox_tons": "NO$_x$ (additive)",
        f"ldn|{COMMUNITY_EDGE}": "$L_{dn}$, Port Isabel edge (energy)",
        f"na100|{COMMUNITY_EDGE}": "Events $>100$ dBA, Port Isabel edge",
        f"days90|{COMMUNITY_EDGE}": "Days $>90$ dBA, Port Isabel edge",
        f"days90|{COMMUNITY_SPI}": "Days $>90$ dBA, S. Padre Island edge",
        f"boom2|{COMMUNITY_EDGE}": "Booms $>2$ psf, Port Isabel edge",
        f"disturbed|{HABITAT}": "Disturbed fraction, habitat (1 h recovery)",
        "closure_h": "Closure hours",
        "weekend_avail": "Weekend daylight availability",
    }
    out = [
        r"\begin{tabular}{lrrrr}",
        r"\hline",
        f"Metric at {int(top)} launches/yr & Simulated & Linear from 25/yr & Error & Relative error (95\\% CI) \\\\",
        r"\hline",
    ]
    for _, r in e.iterrows():
        if r["metric"] not in labels:
            continue
        rel = "" if pd.isna(r["relative_error"]) else f"{100 * r['relative_error']:+.0f}\\%"
        if not pd.isna(r.get("rel_error_lo", float("nan"))):
            rel += f" ({100 * r['rel_error_lo']:+.0f}, {100 * r['rel_error_hi']:+.0f})"
        unit = " dB" if "ldn" in r["metric"] else ""
        out.append(
            f"{labels[r['metric']]} & {fmtnum(float(r['simulated']))} & {fmtnum(float(r['linear_extrapolation']))} & {fmtnum(float(r['error']))}{unit} & {rel} \\\\"
        )
    out += [r"\hline", r"\end{tabular}"]
    return "\n".join(out)


def main() -> None:
    s = summary()
    parts = [
        "% Generated by figures/tables.py from saved results; do not edit by hand.",
        r"\newcommand{\tableScenarios}{" + table_scenarios() + "}",
        r"\newcommand{\tablePerOperation}{" + table_per_operation() + "}",
        r"\newcommand{\tableMetricDefinitions}{" + table_metric_definitions() + "}",
        r"\newcommand{\tableVerification}{" + table_verification() + "}",
        r"\newcommand{\tableAggregationErrors}{" + table_aggregation_errors() + "}",
        f"% pipeline {s.get('pipeline_name')} seed {s.get('seed')}",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(parts) + "\n", encoding="utf-8")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
