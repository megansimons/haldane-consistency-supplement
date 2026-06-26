"""Compute scored real reactions and generate Figures 3 and 4 for the manuscript.

All numbers are taken from primary literature (see manuscript references):
  - PGI kinetics: Stodeman & Schwarz 2004, Anal. Biochem. 329(2), 307-315.
  - PGI thermo:  Goldberg & Tewari 1995 (JPCRD Part 5); TECRDB 60KAH/LOW, 68DYS/NOL.
  - TPI kinetics: Putman, Coulson, Farley, Riddleston & Knowles 1972, Biochem. J. 129, 301-310.
  - TPI thermo:   Veech, Raijman, Dalziel & Krebs 1969, Biochem. J. 115, 837-842 (K=22.0 +/- 0.25).
  - Mandelate racemase: St Maurice & Bearne 2002, Biochemistry 41, 4048-4058.
  - Proline racemase:   Fisher, Albery & Knowles 1986, Biochemistry 25, 2529-2537.
  - Added amino-acid racemases: May et al. 2007 (B. anthracis RacE1/RacE2),
    Gallo & Knowles 1993 (L. fermenti glutamate racemase), Yamashita et al. 2004
    (B. bifidum aspartate racemase), and Strych et al. 2000 (P. aeruginosa DadX).
  - Fumarase thermo:    Goldberg & Tewari 1995 (JPCRD Part 4); Gajewski et al. 1985.
  - Fumarase buffer:    Genda et al. 2006, Biosci. Biotechnol. Biochem. 70, 1102-1109.
"""
from __future__ import annotations

import math
import os
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def J(x: float) -> float:
    return 0.5 * (x + 1.0 / x) - 1.0


def fold(x: float) -> float:
    return x if x >= 1.0 else 1.0 / x


def classify(x: float) -> str:
    f = fold(x)
    if f <= 2:
        return "within twofold"
    if f <= 5:
        return "2-5-fold"
    if f <= 10:
        return "5-10-fold"
    return ">10-fold"


# (label, K_kin, K_thermo, marker-group)
SCORED = [
    ("PGI 60KAH/LOW", 0.307, 0.260, "PGI"),
    ("PGI 68DYS/NOL", 0.395, 0.290, "PGI"),
    ("TPI", 20.4, 22.0, "TPI"),
    ("Mandelate racemase", 1.14, 1.000, "racemase"),
    ("Proline racemase", 1.000, 1.000, "racemase"),
    ("B. anthracis RacE2", (38.0 / 3.7) / (1.6 / 0.20), 1.000, "racemase-added"),
    ("B. anthracis RacE1", (17.7 / 8.0) / (3.9 / 0.90), 1.000, "racemase-added"),
    ("L. fermenti GluR", (69.0 / 0.33) / (68.0 / 0.26), 1.000, "racemase-added"),
    ("B. bifidum AspR", 14.9 / 15.3, 1.000, "racemase-added"),
    ("P. aeruginosa DadX", (155.0 / 1.40) / (134.0 / 1.40), 1.000, "racemase-added"),
]

# Fumarase kinetic Haldane estimates (Genda et al. 2006, C. glutamicum) scored against the
# representative evaluated thermodynamic value K_thermo = 4.2.
# K_kin = (kcat_F K_m,M)/(kcat_M K_m,F).
FUMARASE_KTHERMO = 4.2
FUMARASE = [
    ("Fumarase phosphate pH6", 6.4, FUMARASE_KTHERMO, "fum-P"),
    ("Fumarase phosphate pH7", 6.1, FUMARASE_KTHERMO, "fum-P"),
    ("Fumarase phosphate pH8", 4.6, FUMARASE_KTHERMO, "fum-P"),
    ("Fumarase MES pH6", 16.0, FUMARASE_KTHERMO, "fum-NP"),
    ("Fumarase TES pH7", 19.0, FUMARASE_KTHERMO, "fum-NP"),
    ("Fumarase Tris pH8", 17.0, FUMARASE_KTHERMO, "fum-NP"),
]


def print_table() -> None:
    print(f"{'reaction':22s} {'Kkin':>7s} {'Kthermo':>8s} {'x':>7s} {'fold':>6s} {'C':>8s}  fold band")
    for label, kk, kt, _ in SCORED:
        x = kk / kt
        print(f"{label:22s} {kk:7.3f} {kt:8.3f} {x:7.3f} {fold(x):6.3f} {J(x):8.4f}  {classify(x)}")

    # thermodynamic-side uncertainty from inter-source spread
    pgi_thermo = [0.260, 0.290]
    fum_thermo = [3.98, 4.25, 4.28, 4.43]
    s_pgi = statistics.stdev([math.log(v) for v in pgi_thermo])
    s_fum = statistics.stdev([math.log(v) for v in fum_thermo])
    print(f"\nsigma_ln(thermo) PGI inter-source spread = {s_pgi:.4f}")
    print(f"sigma_ln(thermo) fumarase inter-source spread = {s_fum:.4f}")

    # PGI primary standardized score with both uncertainties
    kk, kt, skk = 0.307, 0.260, 0.053
    lnx = math.log(kk / kt)
    s_kin = skk / kk
    s_delta = math.hypot(s_kin, s_pgi)
    z = lnx / s_delta
    print(f"\nPGI primary: ln x={lnx:.4f}, sigma_kin={s_kin:.4f}, sigma_thermo={s_pgi:.4f},"
          f" sigma_delta={s_delta:.4f}, z={z:.3f}")

    # fumarase kinetic Haldane scores vs thermodynamic K' = 4.2 (Genda 2006)
    print("\nFumarase kinetic Haldane vs thermodynamic K'=4.2 (Genda 2006):")
    for label, kk, kt, _ in FUMARASE:
        x = kk / kt
        print(f"  {label:24s} Kkin={kk:5.1f} x={x:6.3f} fold={fold(x):5.2f}"
              f" C={J(x):8.4f}  {classify(x)}")


def make_consistency_figure(path: str) -> None:
    fig, ax = plt.subplots(figsize=(5.3, 4.6))
    lo, hi = 0.1, 30.0
    grid = [lo, hi]
    ax.plot(grid, grid, color="black", lw=1.2, zorder=1, label="line of identity")
    for fac, c in [(2.0, "0.55"), (5.0, "0.78")]:
        ax.plot([lo, hi], [lo * fac, hi * fac], ls="--", color=c, lw=0.9, zorder=1)
        ax.plot([lo, hi], [lo / fac, hi / fac], ls="--", color=c, lw=0.9, zorder=1)

    groups = {
        "PGI": ("#1f77b4", "o", "PGI (2 source-traced comparators)"),
        "TPI": ("#2ca02c", "s", "TPI"),
        "racemase": ("#d62728", "^", "original racemases ($K'_{thermo}=1$)"),
        "racemase-added": ("#8c564b", "P", "added amino-acid racemases"),
        "fum-P": ("#9467bd", "D", "fumarase (phosphate)"),
        "fum-NP": ("#ff7f0e", "v", "fumarase (non-phosphate)"),
    }
    seen = set()
    for label, kk, kt, grp in SCORED + FUMARASE:
        color, marker, leg = groups[grp]
        ax.scatter([kt], [kk], c=color, marker=marker, s=70, zorder=3,
                   edgecolors="black", linewidths=0.5,
                   label=leg if grp not in seen else None)
        seen.add(grp)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(r"thermodynamic $K'_{\mathrm{eq,thermo}}$")
    ax.set_ylabel(r"kinetic $K'_{\mathrm{eq,kin}}$")
    ax.set_title("Kinetic vs. thermodynamic equilibrium constants")
    ax.text(0.30, 0.55 * 5, "2-, 5-fold bands", fontsize=8, color="0.4", rotation=45)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.95)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def make_distribution_figure(path: str) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.1))
    short = {
        "PGI 60KAH/LOW": "PGI-1",
        "PGI 68DYS/NOL": "PGI-2",
        "TPI": "TPI",
        "Mandelate racemase": "MR",
        "Proline racemase": "PR",
        "B. anthracis RacE2": "RacE2",
        "B. anthracis RacE1": "RacE1",
        "L. fermenti GluR": "Lf-GluR",
        "B. bifidum AspR": "Bb-AspR",
        "P. aeruginosa DadX": "Pa-DadX",
        "Fumarase phosphate pH6": "Fum-P6",
        "Fumarase phosphate pH7": "Fum-P7",
        "Fumarase phosphate pH8": "Fum-P8",
        "Fumarase MES pH6": "Fum-MES6",
        "Fumarase TES pH7": "Fum-TES7",
        "Fumarase Tris pH8": "Fum-Tris8",
    }
    items = [(short[label], J(kk / kt), grp, label) for (label, kk, kt, grp) in SCORED + FUMARASE]
    color_map = {"PGI": "#1f77b4", "TPI": "#2ca02c", "racemase": "#d62728",
                 "racemase-added": "#8c564b", "fum-P": "#9467bd", "fum-NP": "#ff7f0e"}
    labels = [it[0] for it in items]
    floor = 1e-4
    vals = [max(it[1], floor) for it in items]
    colors = [color_map[it[2]] for it in items]

    xs = list(range(len(items)))
    ax.scatter(xs, vals, s=72, color=colors, edgecolor="black", linewidth=0.5, zorder=3)
    ax.set_yscale("log")
    ax.axhline(0.25, ls="--", color="0.4", lw=0.9)
    ax.text(len(items) - 0.5, 0.28, "twofold reporting boundary (0.25)", ha="right",
            va="bottom", fontsize=8, color="0.3")
    for i, (_, cost, _, label) in enumerate(items):
        if cost == 0.0:
            ax.annotate(
                "exact zero\nplotted at floor",
                xy=(i, floor),
                xytext=(i + 0.65, floor * 3.0),
                arrowprops={"arrowstyle": "->", "lw": 0.8, "color": "0.25"},
                fontsize=7,
                color="0.25",
                ha="left",
                va="bottom",
            )
            break
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=7, rotation=35, ha="right")
    ax.set_ylabel(r"Haldane-consistency score $C_{\mathrm{Haldane}}$")
    ax.set_title("Consistency cost for individual scored comparisons")
    ax.set_ylim(floor * 0.7, 2.0)
    ax.grid(axis="y", alpha=0.25, which="both")
    ax.margins(x=0.02)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    print_table()
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
    os.makedirs(out_dir, exist_ok=True)
    make_consistency_figure(os.path.join(out_dir, "fig_consistency.pdf"))
    make_distribution_figure(os.path.join(out_dir, "fig_distribution.pdf"))
    print(f"\nFigures written to {out_dir}")


if __name__ == "__main__":
    main()
