"""Digitize noise contour maps from the 2025 Tiered EA noise appendix.

The KBR noise assessment (Appendix B) presents LAmax, SEL and DNL contours
as colored lines on satellite base maps with a printed scale bar. Each page
is rasterized, the legend colors are matched to their level labels from the
PDF text layer, the scale bar length in pixels is measured, and the pixels
of each contour are converted to polar coordinates about the pad (centroid
of the innermost contour). The radius of each contour is summarized per
azimuth bin and over the overland sector that contains the receptors.

The output is a processed table, not a transcription: it is regenerated
deterministically from the PDF in ``data/raw`` and carries the figure and
page it was derived from.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pymupdf

ROOT = Path(__file__).resolve().parents[3]
APPENDIX_B = ROOT / "data" / "raw" / "faa" / "2025_tiered_ea_appendix_b_noise.pdf"
OUT = ROOT / "data" / "processed" / "digitized_contours.csv"
OUT_BINS = ROOT / "data" / "processed" / "digitized_contours_by_azimuth.csv"

DPI = 300
N_BINS = 16
# Overland sector used for the footprint fits: west through north, where the
# community and habitat receptors lie (azimuth measured clockwise from north).
OVERLAND = (240.0, 360.0)


@dataclass(frozen=True)
class FigureSpec:
    page: int  # 1-based PDF page
    figure: str
    event_type: str
    metric: str
    bar_miles: float
    legend_levels: tuple[float, ...]  # legend swatches in row-major order


LEVELS_SINGLE_EVENT = (90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0)
LEVELS_DNL = (65.0, 85.0, 75.0, 95.0)

FIGURES = [
    FigureSpec(14, "Figure 3", "launch", "lamax", 8.0, LEVELS_SINGLE_EVENT),
    FigureSpec(15, "Figure 4", "launch", "lamax", 4.0, LEVELS_SINGLE_EVENT),
    FigureSpec(16, "Figure 5", "launch", "sel", 8.0, LEVELS_SINGLE_EVENT),
    FigureSpec(17, "Figure 6", "launch", "sel", 4.0, LEVELS_SINGLE_EVENT),
    FigureSpec(21, "Figure 9", "ss_landing", "lamax", 4.0, LEVELS_SINGLE_EVENT),
    FigureSpec(22, "Figure 10", "ss_landing", "sel", 4.0, LEVELS_SINGLE_EVENT),
    FigureSpec(25, "Figure 12", "sh_landing", "lamax", 4.0, LEVELS_SINGLE_EVENT),
    FigureSpec(26, "Figure 13", "sh_landing", "sel", 4.0, LEVELS_SINGLE_EVENT),
    FigureSpec(29, "Figure 15", "static_fire_ss", "lamax", 4.0, LEVELS_SINGLE_EVENT),
    FigureSpec(30, "Figure 16", "static_fire_ss", "sel", 4.0, LEVELS_SINGLE_EVENT),
    FigureSpec(33, "Figure 18", "static_fire", "lamax", 4.0, LEVELS_SINGLE_EVENT),
    FigureSpec(34, "Figure 19", "static_fire", "sel", 4.0, LEVELS_SINGLE_EVENT),
    FigureSpec(37, "Figure 21", "scenario_25_ops", "dnl", 4.0, LEVELS_DNL),
]


def _render(doc: pymupdf.Document, page_no: int) -> tuple[np.ndarray, float]:
    page = doc[page_no - 1]
    zoom = DPI / 72.0
    pix = page.get_pixmap(matrix=pymupdf.Matrix(zoom, zoom), alpha=False)  # type: ignore[no-untyped-call]
    arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)[:, :, :3]
    return arr.astype(np.int16), zoom


def _frame(gray: np.ndarray) -> tuple[int, int, int, int]:
    """Map frame (top, bottom, left, right) from long dark lines."""
    h, w = gray.shape
    dark = gray < 70
    rows = np.flatnonzero(dark.sum(axis=1) > 0.6 * w)
    rows = rows[rows > 0.08 * h]
    top = int(rows[0])
    # the next long line well below the top is the legend separator
    later = rows[rows > top + 0.3 * h]
    bottom = int(later[0]) if later.size else int(0.75 * h)
    cols = np.flatnonzero(dark[top:bottom].sum(axis=0) > 0.5 * (bottom - top))
    left, right = int(cols[0]), int(cols[-1])
    return top, bottom, left, right


def _scale_px_per_mile(gray: np.ndarray, bottom: int, bar_miles: float) -> float:
    """Pixels per mile from the printed scale bar.

    The bar is the longest dark horizontal run in the lower-right part of
    the legend box that is clearly shorter than the box border lines.
    """
    h, w = gray.shape
    x0 = int(0.6 * w)
    region = gray[bottom + 5 : min(h, bottom + int(0.12 * h)), x0 : int(0.9 * w)]
    dark = region < 80
    width = region.shape[1]
    best = 0
    for row in dark:
        run = 0
        longest = 0
        for v in row:
            run = run + 1 if v else 0
            longest = max(longest, run)
        if 0.1 * width < longest < 0.85 * width:
            best = max(best, longest)
    if best < 50:
        raise RuntimeError("scale bar not found")
    return best / bar_miles


def _legend_colors(img: np.ndarray, bottom: int, levels: tuple[float, ...]) -> dict[float, np.ndarray]:
    """Map legend levels to swatch colors.

    The legend is rasterized in the PDF, so the swatches (short saturated
    horizontal line segments in the left part of the legend box) are located
    by image analysis and assigned to ``levels`` in row-major reading order.
    """
    h, w, _ = img.shape
    y0 = bottom + 5
    y1 = min(h, bottom + int(0.12 * h))
    x1 = int(0.45 * w)
    region = img[y0:y1, :x1]
    sat = (region.max(axis=2) - region.min(axis=2)) > 90
    # row profile: rows containing a saturated run of at least 40 px
    swatches = []
    visited = np.zeros_like(sat)
    ys, xs = np.nonzero(sat)
    for y, x in zip(ys, xs, strict=True):
        if visited[y, x]:
            continue
        # flood fill (4-connectivity) within a bounded box
        stack = [(y, x)]
        comp = []
        visited[y, x] = True
        while stack:
            cy, cx = stack.pop()
            comp.append((cy, cx))
            for ny, nx in ((cy + 1, cx), (cy - 1, cx), (cy, cx + 1), (cy, cx - 1)):
                if 0 <= ny < sat.shape[0] and 0 <= nx < sat.shape[1] and sat[ny, nx] and not visited[ny, nx]:
                    visited[ny, nx] = True
                    stack.append((ny, nx))
        cy = np.asarray([c[0] for c in comp])
        cx = np.asarray([c[1] for c in comp])
        if cx.max() - cx.min() >= 40 and cy.max() - cy.min() <= 30:
            color = np.median(region[cy, cx], axis=0)
            swatches.append((float(cy.mean()), float(cx.mean()), color))
    if not swatches:
        return {}
    # group into rows
    swatches.sort(key=lambda t: t[0])
    rows: list[list[tuple[float, float, np.ndarray]]] = []
    for sw in swatches:
        if rows and abs(sw[0] - rows[-1][0][0]) < 25:
            rows[-1].append(sw)
        else:
            rows.append([sw])
    ordered = [sw for row in rows for sw in sorted(row, key=lambda t: t[1])]
    out: dict[float, np.ndarray] = {}
    for lvl, sw in zip(levels, ordered, strict=False):
        out[lvl] = sw[2]
    return out


def _contour_pixels(
    img: np.ndarray, color: np.ndarray, frame: tuple[int, int, int, int], tol: float = 55.0
) -> tuple[np.ndarray, np.ndarray]:
    """Pixels of one contour line: color match, excluding solid fills.

    Contours are thin lines, so a matching pixel whose 15 x 15 neighbourhood
    is mostly matching belongs to a map fill of the same color (for example
    the national seashore polygon) and is discarded.
    """
    from scipy.ndimage import uniform_filter

    top, bottom, left, right = frame
    sub = img[top + 3 : bottom - 3, left + 3 : right - 3]
    d = np.sqrt(((sub - color[None, None, :]) ** 2).sum(axis=2))
    mask = d < tol
    local = uniform_filter(mask.astype(float), size=15, mode="constant")
    mask &= local < 0.45
    ys, xs = np.nonzero(mask)
    return ys + top + 3, xs + left + 3


def digitize_page(doc: pymupdf.Document, spec: FigureSpec) -> tuple[pd.DataFrame, pd.DataFrame]:
    img, zoom = _render(doc, spec.page)
    gray = img.mean(axis=2)
    frame = _frame(gray)
    px_per_mile = _scale_px_per_mile(gray, frame[1], spec.bar_miles)
    colors = _legend_colors(img, frame[1], spec.legend_levels)
    if not colors:
        raise RuntimeError(f"no legend entries found on page {spec.page}")
    pix = {lvl: _contour_pixels(img, c, frame) for lvl, c in colors.items()}
    # pad = centroid of the innermost contour that has enough pixels
    inner = max((lvl for lvl, (ys, xs) in pix.items() if ys.size > 30), default=None)
    if inner is None:
        raise RuntimeError(f"no contour pixels on page {spec.page}")
    ys, xs = pix[inner]
    pad_y, pad_x = float(np.median(ys)), float(np.median(xs))
    bins = np.arange(0, 360 + 1e-9, 360 / N_BINS)
    rows = []
    summary = []
    for lvl in sorted(pix):
        ys, xs = pix[lvl]
        if ys.size < 30:
            continue
        dx = xs - pad_x
        dy = pad_y - ys  # north up
        r = np.hypot(dx, dy) / px_per_mile
        az = np.degrees(np.arctan2(dx, dy)) % 360.0
        keep = r > 0.05
        r, az = r[keep], az[keep]
        # Contours are closed curves of moderate aspect ratio around the pad;
        # pixels far outside the curve (e.g. same-colored map fills) are dropped.
        for _ in range(3):
            med = np.median(r)
            keep = (r < 2.5 * med) & (r > 0.2 * med)
            r, az = r[keep], az[keep]
        per_bin = []
        for b0, b1 in zip(bins[:-1], bins[1:], strict=True):
            m = (az >= b0) & (az < b1)
            if m.sum() >= 5:
                rb = float(np.median(r[m]))
                per_bin.append((0.5 * (b0 + b1), rb, int(m.sum())))
                rows.append(
                    {
                        "page": spec.page,
                        "figure": spec.figure,
                        "event_type": spec.event_type,
                        "metric": spec.metric,
                        "level": lvl,
                        "azimuth_deg": 0.5 * (b0 + b1),
                        "r_miles": rb,
                        "n_pixels": int(m.sum()),
                    }
                )
        if not per_bin:
            continue
        arr = np.asarray(per_bin)
        over = arr[(arr[:, 0] >= OVERLAND[0]) & (arr[:, 0] < OVERLAND[1])]
        summary.append(
            {
                "page": spec.page,
                "figure": spec.figure,
                "event_type": spec.event_type,
                "metric": spec.metric,
                "level": lvl,
                "r_overland_miles": float(np.median(over[:, 1])) if over.size else np.nan,
                "r_overland_min_miles": float(over[:, 1].min()) if over.size else np.nan,
                "r_overland_max_miles": float(over[:, 1].max()) if over.size else np.nan,
                "r_all_min_miles": float(arr[:, 1].min()),
                "r_all_max_miles": float(arr[:, 1].max()),
                "n_bins": int(arr.shape[0]),
                "n_pixels": int(ys.size),
                "px_per_mile": px_per_mile,
                "bar_miles": spec.bar_miles,
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(summary)


def build(pdf: Path = APPENDIX_B) -> pd.DataFrame:
    doc = pymupdf.open(pdf)  # type: ignore[no-untyped-call]
    all_rows, all_sum = [], []
    for spec in FIGURES:
        rows, summ = digitize_page(doc, spec)
        all_rows.append(rows)
        all_sum.append(summ)
    bins_df = pd.concat(all_rows, ignore_index=True)
    summ_df = pd.concat(all_sum, ignore_index=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    bins_df.to_csv(OUT_BINS, index=False, float_format="%.4f")
    summ_df.to_csv(OUT, index=False, float_format="%.4f")
    return summ_df


def as_footprint_points(summary: pd.DataFrame) -> pd.DataFrame:
    """Digitized overland radii in the schema of the transcribed footprint table.

    Where the same level appears on a wide and a zoomed map of the same
    event, the zoomed map (smaller scale bar) is used unless the contour is
    clipped by its frame, in which case the wide map is used.
    """
    rows = []
    s = summary[summary["metric"].isin(["lamax", "sel"]) & summary["r_overland_miles"].notna()].copy()
    for keys, g in s.groupby(["event_type", "metric", "level"]):
        et, met, lvl = str(keys[0]), str(keys[1]), float(str(keys[2]))
        g = g.sort_values("bar_miles")
        pick = None
        for _, r in g.iterrows():
            # a contour that reaches the frame has an overland maximum near the half-width of the map
            if r["r_overland_max_miles"] < 0.95 * r["bar_miles"] * 2.2:
                pick = r
                break
        if pick is None:
            pick = g.iloc[-1]
        pk: dict[str, Any] = {str(k): v for k, v in pick.to_dict().items()}
        rows.append(
            {
                "event_type": et,
                "metric": met,
                "point_id": f"D{int(pk['page'])}_{int(lvl)}",
                "value": float(lvl),
                "value_hi": np.nan,
                "distance_miles": float(pk["r_overland_miles"]),
                "location": np.nan,
                "weight": 1.0,
                "source": f"digitized from 2025 Tiered EA Appendix B {pk['figure']}",
                "page": f"PDF p{int(pk['page'])}",
                "quote": (
                    f"median contour radius over azimuths {OVERLAND[0]:.0f}-{OVERLAND[1]:.0f} deg; "
                    f"{int(pk['n_pixels'])} pixels"
                ),
            }
        )
    return pd.DataFrame(rows)


if __name__ == "__main__":
    print(build().to_string())
