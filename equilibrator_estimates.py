"""eQuilibrator (component-contribution) thermodynamic comparators.

Computes the standard transformed Gibbs energy dG'^o and the corresponding
apparent equilibrium constant K'_eQ for the curated reactions, providing the
*estimated* thermodynamic comparator used in the manuscript's secondary score

    C_Haldane^est = J( K'_kin / K'_eQ ),     J(x) = 1/2 (x + 1/x) - 1.

These are model-based estimates (component contribution), not experimental
ground truth; they are reported alongside, not in place of, the evaluated
experimental comparators (TECRDB).

Requirements
------------
    pip install equilibrator-api

The first invocation downloads the equilibrator compound cache (~hundreds of
MB) via pooch; subsequent runs are fast. Standard conditions used here match
the eQuilibrator defaults reported in the text: pH 7.0, ionic strength
0.25 M, pMg 3.0, T = 298.15 K.
"""
from __future__ import annotations

import math

from equilibrator_api import ComponentContribution, Q_

# eQuilibrator standard conditions (stated in the manuscript).
P_H = 7.0
IONIC_STRENGTH = "0.25M"
P_MG = 3.0
TEMPERATURE_K = 298.15

R_KJ = 8.314462618e-3  # kJ / (mol K)
RT = R_KJ * TEMPERATURE_K  # kJ / mol

# Reactions written in the same direction as the manuscript, with KEGG IDs.
REACTIONS = {
    "PGI  G6P = F6P": "kegg:C00092 = kegg:C00085",
    "TPI  GAP = DHAP": "kegg:C00118 = kegg:C00111",
    "Fumarase  fumarate + H2O = (S)-malate": "kegg:C00122 + kegg:C00001 = kegg:C00149",
}


def _value_error_kJ(measurement) -> tuple[float, float]:
    """Extract (value, error) in kJ/mol from a pint Measurement, robustly."""
    try:
        value = float(measurement.value.m_as("kJ/mol"))
        error = float(measurement.error.m_as("kJ/mol"))
        return value, error
    except AttributeError:
        q = measurement.to("kJ/mol").magnitude
        return float(q.nominal_value), float(q.std_dev)


def main() -> None:
    cc = ComponentContribution()
    cc.p_h = Q_(P_H)
    cc.ionic_strength = Q_(IONIC_STRENGTH)
    cc.p_mg = Q_(P_MG)
    cc.temperature = Q_(f"{TEMPERATURE_K}K")

    print(
        f"eQuilibrator conditions: pH {P_H}, I = {IONIC_STRENGTH}, "
        f"pMg {P_MG}, T = {TEMPERATURE_K} K\n"
    )
    for name, formula in REACTIONS.items():
        rxn = cc.parse_reaction_formula(formula)
        balanced = rxn.is_balanced()
        dg = cc.standard_dg_prime(rxn)
        val, err = _value_error_kJ(dg)
        keq = math.exp(-val / RT)
        flag = "" if balanced else "  [NOT BALANCED]"
        print(
            f"{name:40s} dG'^o = {val:+6.2f} +/- {err:4.2f} kJ/mol   "
            f"K'_eQ = {keq:7.3f}{flag}"
        )


if __name__ == "__main__":
    main()
