"""Figure generation for the Haldane-consistency manuscript.

Produces:

* Figure 1 -- the shape of the reciprocal cost J(x) on a log ratio axis
  (data-free; reproduces the in-text pgfplots figure).
* Figure 3 -- the main consistency scatter, ln K'_eq,kin vs ln K'_eq,thermo,
  with the line of identity and twofold/fivefold/tenfold agreement bands,
  colored by consistency class.
* Figure 4 -- the distribution of C_Haldane across a dataset, optionally
  beside |ln x|.

matplotlib is imported here (and only here); the core ``haldane_consistency``
module stays dependency-light.
"""

from __future__ import annotations

import math
import os
from typing import List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")  # headless / reproducible rendering
import matplotlib.pyplot as plt  # noqa: E402

from haldane_consistency import (  # noqa: E402
    CLASS_LABELS,
    ScoreResult,
    reciprocal_cost,
)

# Consistent palette for the four consistency classes.
_CLASS_COLORS = {
    "consistent": "#2c7bb6",
    "mildly inconsistent": "#abd9e9",
    "strongly inconsistent": "#fdae61",
    "severely inconsistent": "#d7191c",
}

_AXIS_BLUE = "#1f5a9c"


def _logspace(lo: float, hi: float, n: int) -> List[float]:
    """Return ``n`` points geometrically spaced between ``lo`` and ``hi``."""
    if n < 2:
        return [lo]
    ratio = hi / lo
    return [lo * ratio ** (i / (n - 1)) for i in range(n)]


def _ensure_parent_dir(path: str) -> None:
    """Create the parent directory of ``path`` if it does not yet exist."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def figure_cost_shape(
    path: str = "figures/figure1_cost_shape.png",
    x_min: float = 0.05,
    x_max: float = 20.0,
    n: int = 401,
) -> str:
    """Plot J(x) on a log ratio axis and save to ``path``."""
    xs = _logspace(x_min, x_max, n)
    ys = [reciprocal_cost(x) for x in xs]

    fig, ax = plt.subplots(figsize=(7.0, 3.4))
    ax.plot(xs, ys, color=_AXIS_BLUE, lw=2.0)
    ax.set_xscale("log")
    ax.axvline(1.0, ls="--", color="gray", lw=1.0)
    ax.annotate(
        r"$J(10)=4.05$",
        xy=(10.0, reciprocal_cost(10.0)),
        xytext=(10.0, 5.0),
        color=_AXIS_BLUE,
        fontsize=9,
    )
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0.0, 10.0)
    ax.set_xlabel(r"consistency ratio $x = K'_{\mathrm{eq,kin}}/K'_{\mathrm{eq,thermo}}$ (log scale)")
    ax.set_ylabel(r"$J(x) = \cosh(\ln x) - 1$")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    _ensure_parent_dir(path)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def figure_consistency_scatter(
    results: Sequence[ScoreResult],
    path: str = "figures/figure3_consistency_scatter.png",
) -> str:
    """Scatter ln K'_eq,kin vs ln K'_eq,thermo with identity and fold bands."""
    if not results:
        raise ValueError("no results to plot")

    ln_kin = [math.log(r.keq_kin) for r in results]
    ln_thermo = [math.log(r.keq_thermo) for r in results]
    classes = [r.consistency_class for r in results]

    lo = min(min(ln_kin), min(ln_thermo))
    hi = max(max(ln_kin), max(ln_thermo))
    pad = 0.1 * (hi - lo) if hi > lo else 1.0
    lo, hi = lo - pad, hi + pad

    fig, ax = plt.subplots(figsize=(5.2, 5.0))

    # Identity line and symmetric agreement bands (twofold, fivefold, tenfold).
    ax.plot([lo, hi], [lo, hi], color="black", lw=1.2, label="identity")
    band_styles = [
        (math.log(2.0), "#2c7bb6", "twofold"),
        (math.log(5.0), "#fdae61", "fivefold"),
        (math.log(10.0), "#d7191c", "tenfold"),
    ]
    for width, color, label in band_styles:
        ax.fill_between(
            [lo, hi],
            [lo - width, hi - width],
            [lo + width, hi + width],
            color=color,
            alpha=0.10,
            label=f"within {label}",
        )

    # Points colored by class.
    for cls in CLASS_LABELS:
        xs = [t for t, c in zip(ln_thermo, classes) if c == cls]
        ys = [k for k, c in zip(ln_kin, classes) if c == cls]
        if xs:
            ax.scatter(
                xs,
                ys,
                s=42,
                color=_CLASS_COLORS[cls],
                edgecolor="black",
                linewidth=0.4,
                zorder=3,
                label=cls,
            )

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$\ln K'_{\mathrm{eq,thermo}}$")
    ax.set_ylabel(r"$\ln K'_{\mathrm{eq,kin}}$")
    # Place the legend in the lower-right, away from the upper-left region where
    # over-estimating (large C) records cluster.
    ax.legend(fontsize=7, loc="lower right", framealpha=0.9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    _ensure_parent_dir(path)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def figure_score_distribution(
    results: Sequence[ScoreResult],
    path: str = "figures/figure4_score_distribution.png",
    with_log_error: bool = True,
) -> str:
    """Histogram of C_Haldane, optionally beside |ln x|."""
    if not results:
        raise ValueError("no results to plot")

    cost = [r.C_haldane for r in results]

    if with_log_error:
        log_err = [r.abs_log_error for r in results]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.0, 3.2))
        ax1.hist(cost, bins=12, color=_AXIS_BLUE, edgecolor="white")
        ax1.set_xlabel(r"$C_{\mathrm{Haldane}} = \cosh(\ln x) - 1$")
        ax1.set_ylabel("count")
        ax2.hist(log_err, bins=12, color="#7b9fc4", edgecolor="white")
        ax2.set_xlabel(r"$|\ln x|$")
        ax2.set_ylabel("count")
    else:
        fig, ax1 = plt.subplots(figsize=(5.0, 3.2))
        ax1.hist(cost, bins=12, color=_AXIS_BLUE, edgecolor="white")
        ax1.set_xlabel(r"$C_{\mathrm{Haldane}} = \cosh(\ln x) - 1$")
        ax1.set_ylabel("count")

    fig.tight_layout()
    _ensure_parent_dir(path)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


if __name__ == "__main__":
    # Running this module directly regenerates the data-free Figure 1.
    out = figure_cost_shape()
    print(f"wrote {out}")
