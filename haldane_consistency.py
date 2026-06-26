"""Reciprocal Haldane-consistency cost for reversible enzyme kinetics.

Reference implementation of the scoring described in

    "A Reciprocal Reporting Scale for Haldane Consistency in Reversible
     Enzyme Kinetics: A Curated Proof of Concept."

The central quantity is the reciprocal recognition cost

    J(x) = 1/2 (x + 1/x) - 1 = cosh(ln x) - 1,        (Section 3.3)

evaluated at the ratio of the kinetic and thermodynamic estimates of the
apparent (transformed) equilibrium constant,

    C_Haldane = J( K'_eq,kin / K'_eq,thermo ).

This module is intentionally dependency-light: it uses only the Python
standard library (``math``, ``csv``, ``dataclasses``). Plotting lives in the
companion ``figures.py`` module, which is the only place that needs
matplotlib.

Equation/section labels in the docstrings refer to the manuscript.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, fields
from typing import Iterable, List, Optional

__all__ = [
    "R_GAS",
    "T_DEFAULT",
    "reciprocal_cost",
    "reciprocal_cost_from_log",
    "keq_kinetic",
    "ln_keq_kinetic",
    "keq_thermo_from_dg",
    "dg_from_keq",
    "consistency_ratio",
    "haldane_score",
    "fold_error",
    "classify",
    "sigma_ln_keq_kin",
    "sigma_delta",
    "standardized_residual",
    "baseline_scores",
    "EnzymeRecord",
    "ScoreResult",
    "score_record",
    "read_records",
    "write_results",
]

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

R_GAS: float = 8.314462618   # molar gas constant, J / (mol K)
T_DEFAULT: float = 298.15    # default temperature, K (25 degrees C)

# Fixed fold-error reporting bands. These are descriptive labels, not
# validated consistency classes or decision thresholds.
CLASS_BOUNDS = (2.0, 5.0, 10.0)
CLASS_LABELS = (
    "within twofold",
    "2-5-fold",
    "5-10-fold",
    ">10-fold",
)


# ---------------------------------------------------------------------------
# Core cost function (Section 3.3)
# ---------------------------------------------------------------------------

def reciprocal_cost(x: float) -> float:
    """Return the reciprocal recognition cost ``J(x) = 1/2 (x + 1/x) - 1``.

    Defined for ``x > 0``. It is non-negative, zero only at ``x = 1``,
    symmetric under ``x -> 1/x``, and equals ``cosh(ln x) - 1``.
    """
    if x <= 0:
        raise ValueError(f"consistency ratio must be positive, got {x!r}")
    return 0.5 * (x + 1.0 / x) - 1.0


def reciprocal_cost_from_log(delta: float) -> float:
    """Return the cost in the log-ratio coordinate, ``cosh(delta) - 1``.

    Here ``delta = ln x``. Numerically preferable when ``x`` spans many orders
    of magnitude.
    """
    return math.cosh(delta) - 1.0


# ---------------------------------------------------------------------------
# Haldane relation and the thermodynamic estimator (Sections 3.1, 3.2)
# ---------------------------------------------------------------------------

def keq_kinetic(kcat_f: float, km_s: float, kcat_r: float, km_p: float) -> float:
    """Kinetic equilibrium constant from the uni--uni Haldane relation.

    ``K'_eq,kin = (kcat+ * K_M,P) / (kcat- * K_M,S)``  (Section 3.1).

    The kinetic constants must be positive and expressed in consistent units;
    the turnover units cancel, as do the Michaelis-constant units, so the
    result is dimensionless for an ``S <-> P`` reaction.
    """
    _require_positive(kcat_f=kcat_f, km_s=km_s, kcat_r=kcat_r, km_p=km_p)
    return (kcat_f * km_p) / (kcat_r * km_s)


def ln_keq_kinetic(kcat_f: float, km_s: float, kcat_r: float, km_p: float) -> float:
    """Log kinetic equilibrium constant (Section 3.5, additive decomposition).

    ``ln K'_eq,kin = ln kcat+ + ln K_M,P - ln kcat- - ln K_M,S``.
    """
    _require_positive(kcat_f=kcat_f, km_s=km_s, kcat_r=kcat_r, km_p=km_p)
    return math.log(kcat_f) + math.log(km_p) - math.log(kcat_r) - math.log(km_s)


def keq_thermo_from_dg(dg_J_per_mol: float, temperature_K: float = T_DEFAULT) -> float:
    """Thermodynamic equilibrium constant from a standard transformed Gibbs energy.

    ``K'_eq,thermo = exp(-dG'/(R T))``  (Section 3.2). ``dg_J_per_mol`` is in
    joules per mole.
    """
    return math.exp(-dg_J_per_mol / (R_GAS * temperature_K))


def dg_from_keq(keq: float, temperature_K: float = T_DEFAULT) -> float:
    """Standard transformed Gibbs energy from an equilibrium constant.

    ``dG' = -R T ln K'`` in joules per mole.
    """
    if keq <= 0:
        raise ValueError(f"equilibrium constant must be positive, got {keq!r}")
    return -R_GAS * temperature_K * math.log(keq)


# ---------------------------------------------------------------------------
# Score and classification
# ---------------------------------------------------------------------------

def consistency_ratio(keq_kin: float, keq_thermo: float) -> float:
    """Return the consistency ratio ``x = K'_eq,kin / K'_eq,thermo``."""
    if keq_thermo <= 0:
        raise ValueError(f"thermodynamic K' must be positive, got {keq_thermo!r}")
    return keq_kin / keq_thermo


def haldane_score(keq_kin: float, keq_thermo: float) -> float:
    """Return ``C_Haldane = J(K'_eq,kin / K'_eq,thermo)``."""
    return reciprocal_cost(consistency_ratio(keq_kin, keq_thermo))


def fold_error(x: float) -> float:
    """Symmetric fold error ``max(x, 1/x)``; always >= 1."""
    if x <= 0:
        raise ValueError(f"ratio must be positive, got {x!r}")
    return max(x, 1.0 / x)


def classify(x: float) -> str:
    """Assign the fixed fold-error reporting band from the ratio ``x``."""
    f = fold_error(x)
    if f <= CLASS_BOUNDS[0]:
        return CLASS_LABELS[0]
    if f <= CLASS_BOUNDS[1]:
        return CLASS_LABELS[1]
    if f <= CLASS_BOUNDS[2]:
        return CLASS_LABELS[2]
    return CLASS_LABELS[3]


# ---------------------------------------------------------------------------
# Uncertainty-aware score (Section 3.5)
# ---------------------------------------------------------------------------

def sigma_ln_keq_kin(
    sigma_ln_kcat_f: float,
    sigma_ln_km_p: float,
    sigma_ln_kcat_r: float,
    sigma_ln_km_s: float,
) -> float:
    """Standard deviation of ``ln K'_eq,kin`` under independent errors.

    ``sigma^2 = sigma^2_{ln kcat+} + sigma^2_{ln K_M,P}
               + sigma^2_{ln kcat-} + sigma^2_{ln K_M,S}``  (Section 3.5).
    """
    return math.sqrt(
        sigma_ln_kcat_f ** 2
        + sigma_ln_km_p ** 2
        + sigma_ln_kcat_r ** 2
        + sigma_ln_km_s ** 2
    )


def sigma_delta(sigma_ln_kin: float, sigma_ln_thermo: float) -> float:
    """Combined log-ratio standard deviation ``sqrt(sigma^2_kin + sigma^2_thermo)``."""
    return math.hypot(sigma_ln_kin, sigma_ln_thermo)


def standardized_residual(delta: float, s_delta: float) -> float:
    """Standardized residual ``z = delta / sigma_delta`` (Section 3.5)."""
    if s_delta <= 0:
        raise ValueError(f"sigma_delta must be positive, got {s_delta!r}")
    return delta / s_delta


# ---------------------------------------------------------------------------
# Baseline discrepancy measures (Section 4.7)
# ---------------------------------------------------------------------------

def baseline_scores(x: float, temperature_K: float = T_DEFAULT) -> dict:
    """Return the common alternative discrepancy measures for ratio ``x``.

    Includes the absolute fold error ``|x - 1|``, the log error ``|ln x|``,
    the squared log error ``(ln x)^2``, the free-energy error
    ``|ddG| = R T |ln x|`` (in kJ/mol), its square, and the reciprocal cost
    itself, for side-by-side comparison (manuscript scoring-function comparison table).
    """
    delta = math.log(x)
    rt = R_GAS * temperature_K
    ddg_kJ = rt * abs(delta) / 1000.0
    return {
        "abs_fold_error": abs(x - 1.0),
        "log_error": abs(delta),
        "squared_log_error": delta ** 2,
        "free_energy_error_kJ_per_mol": ddg_kJ,
        "squared_free_energy_error_kJ2": ddg_kJ ** 2,
        "haldane_score": reciprocal_cost(x),
    }


# ---------------------------------------------------------------------------
# Record-level scoring
# ---------------------------------------------------------------------------

@dataclass
class EnzymeRecord:
    """A single reversible uni--uni record (one row of Dataset A).

    Required for scoring: ``kcat_f``, ``km_s``, ``kcat_r``, ``km_p`` and
    exactly one thermodynamic comparator (``keq_thermo`` or
    ``dg_thermo_kJ_per_mol``). All ``sigma_ln_*`` fields are optional
    1-sigma uncertainties in natural-log units; if all of them are present
    the standardized score is also computed.
    """

    record_id: str = ""
    ec_number: str = ""
    enzyme_name: str = ""
    organism: str = ""
    reaction_rhea: str = ""
    kinetic_source: str = ""
    thermo_source: str = ""
    pH: Optional[float] = None
    temperature_K: Optional[float] = None
    # kinetic constants: kcat in 1/s, K_M in mM (units cancel in the ratio)
    kcat_f: Optional[float] = None
    km_s: Optional[float] = None
    kcat_r: Optional[float] = None
    km_p: Optional[float] = None
    # thermodynamic comparator: provide one of the two
    keq_thermo: Optional[float] = None
    dg_thermo_kJ_per_mol: Optional[float] = None
    # optional log-scale 1-sigma uncertainties
    sigma_ln_kcat_f: Optional[float] = None
    sigma_ln_km_s: Optional[float] = None
    sigma_ln_kcat_r: Optional[float] = None
    sigma_ln_km_p: Optional[float] = None
    sigma_ln_keq_thermo: Optional[float] = None


@dataclass
class ScoreResult:
    """Computed quantities for one record."""

    record_id: str
    keq_kin: float
    keq_thermo: float
    ratio_x: float
    ln_ratio: float
    C_haldane: float
    consistency_class: str
    z: Optional[float] = None
    # baseline measures
    abs_log_error: Optional[float] = None
    squared_log_error: Optional[float] = None
    abs_ddg_kJ_per_mol: Optional[float] = None


def score_record(rec: EnzymeRecord) -> ScoreResult:
    """Score a single :class:`EnzymeRecord`.

    Raises ``ValueError`` if the kinetic constants are missing/non-positive or
    if no thermodynamic comparator is supplied.
    """
    if None in (rec.kcat_f, rec.km_s, rec.kcat_r, rec.km_p):
        raise ValueError(
            f"record {rec.record_id!r} is missing one or more kinetic constants"
        )

    temperature = rec.temperature_K if rec.temperature_K else T_DEFAULT

    keq_kin = keq_kinetic(rec.kcat_f, rec.km_s, rec.kcat_r, rec.km_p)

    if rec.keq_thermo is not None:
        keq_thermo = rec.keq_thermo
    elif rec.dg_thermo_kJ_per_mol is not None:
        keq_thermo = keq_thermo_from_dg(rec.dg_thermo_kJ_per_mol * 1000.0, temperature)
    else:
        raise ValueError(
            f"record {rec.record_id!r} needs keq_thermo or dg_thermo_kJ_per_mol"
        )

    x = consistency_ratio(keq_kin, keq_thermo)
    delta = math.log(x)
    cost = reciprocal_cost_from_log(delta)
    cls = classify(x)

    # Optional standardized residual, only when all uncertainties are present.
    z = None
    kin_sigmas = (
        rec.sigma_ln_kcat_f,
        rec.sigma_ln_km_p,
        rec.sigma_ln_kcat_r,
        rec.sigma_ln_km_s,
    )
    if all(s is not None for s in kin_sigmas) and rec.sigma_ln_keq_thermo is not None:
        s_kin = sigma_ln_keq_kin(*kin_sigmas)
        s_delta = sigma_delta(s_kin, rec.sigma_ln_keq_thermo)
        if s_delta > 0:
            z = standardized_residual(delta, s_delta)

    rt_kJ = R_GAS * temperature / 1000.0
    return ScoreResult(
        record_id=rec.record_id,
        keq_kin=keq_kin,
        keq_thermo=keq_thermo,
        ratio_x=x,
        ln_ratio=delta,
        C_haldane=cost,
        consistency_class=cls,
        z=z,
        abs_log_error=abs(delta),
        squared_log_error=delta ** 2,
        abs_ddg_kJ_per_mol=rt_kJ * abs(delta),
    )


def score_records(records: Iterable[EnzymeRecord]) -> List[ScoreResult]:
    """Score an iterable of records, returning a list of results."""
    return [score_record(rec) for rec in records]


# ---------------------------------------------------------------------------
# CSV input / output (standard library only)
# ---------------------------------------------------------------------------

# Fields parsed as floats when reading the input CSV.
_FLOAT_FIELDS = {
    "pH",
    "temperature_K",
    "kcat_f",
    "km_s",
    "kcat_r",
    "km_p",
    "keq_thermo",
    "dg_thermo_kJ_per_mol",
    "sigma_ln_kcat_f",
    "sigma_ln_km_s",
    "sigma_ln_kcat_r",
    "sigma_ln_km_p",
    "sigma_ln_keq_thermo",
}

_RECORD_FIELD_NAMES = {f.name for f in fields(EnzymeRecord)}


def _parse_float(value: Optional[str]) -> Optional[float]:
    """Parse a CSV cell into a float, treating blanks as ``None``."""
    if value is None:
        return None
    value = value.strip()
    if value == "" or value.upper() in {"NA", "NAN", "NULL", "NONE"}:
        return None
    return float(value)


def _require_positive(**values: float) -> None:
    for name, val in values.items():
        if val is None or val <= 0:
            raise ValueError(f"{name} must be positive, got {val!r}")


def read_records(path: str) -> List[EnzymeRecord]:
    """Read records from a CSV file with a header row.

    Unknown columns are ignored; missing optional columns default to ``None``
    or the empty string. Numeric columns are coerced to floats.
    """
    records: List[EnzymeRecord] = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            kwargs = {}
            for key, value in raw.items():
                if key is None or key not in _RECORD_FIELD_NAMES:
                    continue
                if key in _FLOAT_FIELDS:
                    kwargs[key] = _parse_float(value)
                else:
                    kwargs[key] = (value or "").strip()
            records.append(EnzymeRecord(**kwargs))
    return records


# Columns written to the scored output CSV (input identity + computed fields).
_OUTPUT_COLUMNS = [
    "record_id",
    "ec_number",
    "enzyme_name",
    "organism",
    "reaction_rhea",
    "kinetic_source",
    "thermo_source",
    "pH",
    "temperature_K",
    "kcat_f",
    "km_s",
    "kcat_r",
    "km_p",
    "keq_kin",
    "keq_thermo",
    "ratio_x",
    "ln_ratio",
    "C_haldane",
    "z",
    "consistency_class",
    "abs_log_error",
    "squared_log_error",
    "abs_ddg_kJ_per_mol",
]


def _fmt(value: Optional[float], digits: int = 6) -> str:
    """Format a float for CSV output; blank for ``None``."""
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def write_results(
    path: str,
    records: Iterable[EnzymeRecord],
    results: Iterable[ScoreResult],
) -> None:
    """Write a CSV that merges the input identity columns with the scores."""
    rec_by_id = {r.record_id: r for r in records}
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(_OUTPUT_COLUMNS)
        for res in results:
            rec = rec_by_id.get(res.record_id, EnzymeRecord())
            writer.writerow(
                [
                    res.record_id,
                    rec.ec_number,
                    rec.enzyme_name,
                    rec.organism,
                    rec.reaction_rhea,
                    rec.kinetic_source,
                    rec.thermo_source,
                    _fmt(rec.pH, 4),
                    _fmt(rec.temperature_K, 6),
                    _fmt(rec.kcat_f),
                    _fmt(rec.km_s),
                    _fmt(rec.kcat_r),
                    _fmt(rec.km_p),
                    _fmt(res.keq_kin),
                    _fmt(res.keq_thermo),
                    _fmt(res.ratio_x),
                    _fmt(res.ln_ratio),
                    _fmt(res.C_haldane),
                    _fmt(res.z),
                    res.consistency_class,
                    _fmt(res.abs_log_error),
                    _fmt(res.squared_log_error),
                    _fmt(res.abs_ddg_kJ_per_mol),
                ]
            )
