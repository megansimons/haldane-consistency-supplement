"""Tests for the Haldane-consistency reference implementation.

Verifies the mathematical properties claimed in the manuscript and reproduces
the values in the worked illustrative example. Runs under pytest, or
standalone via ``python test_haldane_consistency.py``.
"""

from __future__ import annotations

import math

import haldane_consistency as hc

TOL = 1e-9


def _close(a, b, tol=1e-6):
    return abs(a - b) <= tol


# --- core cost properties (Section 3.3) -----------------------------------

def test_normalization_and_nonnegativity():
    assert hc.reciprocal_cost(1.0) == 0.0
    for x in (0.01, 0.3, 1.0, 2.0, 7.5, 100.0):
        assert hc.reciprocal_cost(x) >= -TOL


def test_reciprocal_symmetry():
    for x in (1.3, 2.0, 5.0, 10.0, 42.0):
        assert _close(hc.reciprocal_cost(x), hc.reciprocal_cost(1.0 / x))


def test_tabulated_values():
    # Property 5 / the manuscript's fold-error interpretation of J.
    assert _close(hc.reciprocal_cost(2.0), 0.25)
    assert _close(hc.reciprocal_cost(3.0), 2.0 / 3.0, tol=1e-3)
    assert _close(hc.reciprocal_cost(5.0), 1.6)
    assert _close(hc.reciprocal_cost(10.0), 4.05)
    assert _close(hc.reciprocal_cost(100.0), 49.005)


def test_cosh_form_matches():
    for x in (0.2, 0.9, 1.0, 3.3, 25.0):
        assert _close(hc.reciprocal_cost(x), hc.reciprocal_cost_from_log(math.log(x)))
        assert _close(hc.reciprocal_cost(x), math.cosh(math.log(x)) - 1.0)


def test_small_error_quadratic():
    # C = delta^2/2 + delta^4/24 + O(delta^6) for small delta (Property 4).
    for delta in (1e-3, 5e-3, 1e-2, 0.05):
        quadratic = 0.5 * delta ** 2
        quartic = delta ** 4 / 24.0
        exact = hc.reciprocal_cost_from_log(delta)
        # residual beyond the quartic term is O(delta^6); check it is negligible.
        assert abs(exact - quadratic - quartic) < 1e-6 * quadratic


def test_negative_ratio_raises():
    for bad in (0.0, -1.0):
        try:
            hc.reciprocal_cost(bad)
        except ValueError:
            pass
        else:  # pragma: no cover
            raise AssertionError("expected ValueError for non-positive ratio")


# --- Haldane relation and thermodynamics (Sections 3.1, 3.2) --------------

def test_haldane_relation():
    # R2 of the illustrative example: (200 * 0.4) / (100 * 0.5) = 1.6
    assert _close(hc.keq_kinetic(200, 0.5, 100, 0.4), 1.6)


def test_ln_keq_matches_log_of_keq():
    keq = hc.keq_kinetic(300, 0.1, 30, 0.25)
    assert _close(math.log(keq), hc.ln_keq_kinetic(300, 0.1, 30, 0.25))


def test_dg_keq_roundtrip():
    for keq in (0.01, 1.0, 55.0):
        dg = hc.dg_from_keq(keq, hc.T_DEFAULT)
        assert _close(hc.keq_thermo_from_dg(dg, hc.T_DEFAULT), keq, tol=1e-6 * keq + 1e-9)


# --- fold-band reporting labels ------------------------------------------

def test_classification_boundaries():
    assert hc.classify(1.0) == "within twofold"
    assert hc.classify(2.0) == "within twofold"
    assert hc.classify(0.5) == "within twofold"
    assert hc.classify(3.0) == "2-5-fold"
    assert hc.classify(1.0 / 3.0) == "2-5-fold"
    assert hc.classify(5.0) == "2-5-fold"
    assert hc.classify(8.0) == "5-10-fold"
    assert hc.classify(10.0) == "5-10-fold"
    assert hc.classify(25.0) == ">10-fold"
    assert hc.classify(0.02) == ">10-fold"


def test_fold_error_symmetry():
    assert _close(hc.fold_error(4.0), hc.fold_error(0.25))
    assert hc.fold_error(1.0) == 1.0


# --- uncertainty propagation (Section 3.5) --------------------------------

def test_sigma_propagation():
    s = hc.sigma_ln_keq_kin(0.3, 0.3, 0.3, 0.3)
    assert _close(s, math.sqrt(4 * 0.09))  # 0.6
    sd = hc.sigma_delta(0.6, 0.2)
    assert _close(sd, math.hypot(0.6, 0.2))


def test_standardized_residual():
    z = hc.standardized_residual(1.2, 0.6)
    assert _close(z, 2.0)


# --- baselines (Section 4.7) ----------------------------------------------

def test_baselines():
    b = hc.baseline_scores(10.0, hc.T_DEFAULT)
    assert _close(b["log_error"], math.log(10.0))
    assert _close(b["squared_log_error"], math.log(10.0) ** 2)
    # |ddG| = R T |ln x| in kJ/mol
    assert _close(b["free_energy_error_kJ_per_mol"],
                  hc.R_GAS * hc.T_DEFAULT * math.log(10.0) / 1000.0)
    assert _close(b["haldane_score"], 4.05)


# --- worked illustrative example ------------------------------------------

EXPECTED = {
    # record_id: (keq_kin, x, C_haldane, class)
    "R1": (1.00, 1.0000, 0.00000, "within twofold"),
    "R2": (1.60, 1.0667, 0.00208, "within twofold"),
    "R3": (4.50, 3.0000, 0.66667, "2-5-fold"),
    "R4": (25.0, 8.3333, 3.22667, "5-10-fold"),
    "R5": (50.0, 25.000, 11.5200, ">10-fold"),
}


def _example_records():
    return [
        hc.EnzymeRecord(record_id="R1", kcat_f=100, km_s=0.50, kcat_r=80, km_p=0.40,
                        keq_thermo=1.00, temperature_K=298.15),
        hc.EnzymeRecord(record_id="R2", kcat_f=200, km_s=0.50, kcat_r=100, km_p=0.40,
                        keq_thermo=1.50, temperature_K=298.15),
        hc.EnzymeRecord(record_id="R3", kcat_f=150, km_s=0.20, kcat_r=50, km_p=0.30,
                        keq_thermo=1.50, temperature_K=298.15),
        hc.EnzymeRecord(record_id="R4", kcat_f=300, km_s=0.10, kcat_r=30, km_p=0.25,
                        keq_thermo=3.00, temperature_K=298.15),
        hc.EnzymeRecord(record_id="R5", kcat_f=500, km_s=0.05, kcat_r=20, km_p=0.10,
                        keq_thermo=2.00, temperature_K=298.15),
    ]


def test_illustrative_example_values():
    for rec in _example_records():
        res = hc.score_record(rec)
        keq, x, cost, cls = EXPECTED[rec.record_id]
        assert _close(res.keq_kin, keq, tol=1e-2)
        assert _close(res.ratio_x, x, tol=1e-3)
        assert _close(res.C_haldane, cost, tol=1e-3)
        assert res.consistency_class == cls


def test_score_record_requires_thermo():
    rec = hc.EnzymeRecord(record_id="bad", kcat_f=1, km_s=1, kcat_r=1, km_p=1)
    try:
        hc.score_record(rec)
    except ValueError:
        pass
    else:  # pragma: no cover
        raise AssertionError("expected ValueError when no thermo comparator given")


def test_standardized_residual_only_with_full_uncertainty():
    full = hc.EnzymeRecord(record_id="u", kcat_f=200, km_s=0.5, kcat_r=100, km_p=0.4,
                           keq_thermo=1.5, sigma_ln_kcat_f=0.3, sigma_ln_km_s=0.3,
                           sigma_ln_kcat_r=0.3, sigma_ln_km_p=0.3, sigma_ln_keq_thermo=0.2)
    res_full = hc.score_record(full)
    assert res_full.z is not None

    partial = hc.EnzymeRecord(record_id="p", kcat_f=200, km_s=0.5, kcat_r=100, km_p=0.4,
                              keq_thermo=1.5)
    res_partial = hc.score_record(partial)
    assert res_partial.z is None


# --- standalone runner -----------------------------------------------------

def _run_all() -> int:
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    failed = 0
    for test in tests:
        try:
            test()
            print(f"PASS  {test.__name__}")
        except AssertionError as exc:
            failed += 1
            print(f"FAIL  {test.__name__}: {exc}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"ERROR {test.__name__}: {exc!r}")
    print(f"\n{len(tests) - failed}/{len(tests)} tests passed.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_run_all())
