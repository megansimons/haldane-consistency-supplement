# Haldane Consistency as a Reciprocal Reporting Scale for Reversible Enzyme Kinetics: A Curated Benchmark

Reference code and supplemental data for the manuscript **"Haldane Consistency
as a Reciprocal Reporting Scale for Reversible Enzyme Kinetics: A Curated
Benchmark."**

Archived release: the concept DOI
[10.5281/zenodo.20790110](https://doi.org/10.5281/zenodo.20790110) always resolves to the latest
archived version (currently `v1.4`); each tagged release also receives its own immutable version
DOI on Zenodo.
The live repository is
[megansimons/haldane-consistency-supplement](https://github.com/megansimons/haldane-consistency-supplement).

For a reversible uni--uni reaction `S <-> P`, the forward and reverse kinetic
constants fix an apparent equilibrium constant through the Haldane relation,
which should agree with the value implied by biochemical thermodynamics. The
score quantifies their disagreement with the reciprocal recognition cost

```
C_Haldane = J( K'_eq,kin / K'_eq,thermo ),    J(x) = 1/2 (x + 1/x) - 1 = cosh(ln x) - 1
```

where

```
K'_eq,kin    = (kcat+ * K_M,P) / (kcat- * K_M,S)        (Haldane relation, Section 3.1)
K'_eq,thermo = exp(-dG'/(R T))   or a measured K'        (Section 3.2)
```

`J` is zero exactly at agreement, symmetric under `x -> 1/x`, and equal to a
`cosh` of the free-energy discrepancy in units of `RT`.

## Contents

| File | Purpose |
|------|---------|
| `haldane_consistency.py` | Core library: the cost `J`, the Haldane relation, the thermodynamic estimator, uncertainty propagation, fold-band labeling, baselines, and CSV I/O. **Standard library only.** |
| `curated_reactions.csv` | **The curated real-reaction dataset.** One row per kinetic--thermodynamic comparison, with the kinetic estimate `keq_kin`, the thermodynamic comparator `keq_thermo`, conditions, and primary-literature provenance. |
| `curated_reactions_scored.csv` | Scored output from `score_curated.py`, including `x`, `C_haldane`, fold error, and the fixed fold-discrepancy band. |
| `strenda_metadata.csv` | STRENDA-aligned conditions and directional kinetic terms for the records summarized in the manuscript's STRENDA tables. Missing fields are retained as `not reported`. |
| `Supplementary_Table_S1_reaction_normalization.csv` | Explicit Supplementary Table S1 export: normalized reaction identities, reaction orientation, thermodynamic anchor/comparator, and Rhea/ChEBI audit status. |
| `Supplementary_Table_S2_kinetic_metadata.csv` | Explicit Supplementary Table S2 export: STRENDA-style kinetic metadata for directly documented records. |
| `Supplementary_Table_S3_provenance_limited_controls.csv` | Explicit Supplementary Table S3 export: provenance-limited racemase controls and their interpretation limits. |
| `score_curated.py` | Scores `curated_reactions.csv`, prints the summary table, and writes a scored CSV. **Standard library only.** |
| `equilibrator_estimates.py` | Recomputes the eQuilibrator (component-contribution) thermodynamic comparators used for the secondary score `C_Haldane^est`. Requires `equilibrator-api` (optional). |
| `make_real_figures.py` | Generates the manuscript's real-data Figure 3 (consistency scatter) and Figure 4 (score distribution) into `figures/`, and prints the per-reaction scores and uncertainty analysis. Requires `matplotlib`. |
| `score_curated_output.txt`, `make_real_figures_output.txt` | Captured stdout from the regenerated scoring and figure-output commands. |
| `EXPECTED_OUTPUT_HASHES.sha256` | SHA-256 manifest for regenerated outputs. |
| `RELEASE_VERIFICATION.md` | v1.4 local verification record and post-upload DOI/tag checklist. |
| `figures.py` | The data-free Figure 1 (cost shape) and schematic figure helpers. Requires `matplotlib`. |
| `run_analysis.py` | Command-line driver for the four-constant input format: score a CSV, write an augmented CSV, print a summary, optionally make figures. |
| `example_dataset.csv` | The synthetic worked example (records R1--R5 of the manuscript's worked illustrative example). **Illustrative only; not experimental data.** |
| `test_haldane_consistency.py` | Tests for the mathematical properties and the worked illustrative example values. |
| `lean_formalization_details.md` | Pinned repository, module, and declaration names for an optional, independent Lean 4 formalization of the d'Alembert/cosh classification. The manuscript proof is self-contained (Appendix A, Lemma 3); this material is **not required by or cited in the manuscript** and is provided only as a supplementary independent check. |
| `requirements.txt` | `matplotlib` (figures); `equilibrator-api` (optional, for `equilibrator_estimates.py`). |

## Requirements

* Python 3.8+ (core library, CLI, and tests need nothing beyond the standard library).
* `matplotlib >= 3.3` for figures: `pip install -r requirements.txt`.

## Quick start

Score the bundled synthetic example and write the results:

```bash
python run_analysis.py --input example_dataset.csv --output scored_example.csv
```

Also (re)generate the figures into `figures/`:

```bash
python run_analysis.py -i example_dataset.csv -o scored_example.csv --figures
```

Regenerate the data-free Figure 1 on its own:

```bash
python figures.py
```

Run the tests:

```bash
python test_haldane_consistency.py      # standalone
# or, if pytest is installed:
pytest -q
```

## Reproducing the manuscript's real-data results

Score the curated real-reaction dataset and write the scored CSV:

```bash
python score_curated.py --output curated_reactions_scored.csv
```

Regenerate the real-data Figures 3 and 4 (and print the per-reaction scores and
uncertainty analysis) into `figures/`:

```bash
python make_real_figures.py
```

Recompute the eQuilibrator (component-contribution) thermodynamic comparators
for the secondary score `C_Haldane^est` (optional; needs `equilibrator-api` and
downloads a ~1.3 GB compound cache on first use):

```bash
pip install equilibrator-api
python equilibrator_estimates.py
```

The curated dataset is the authoritative record of where each number comes
from. Each row of `curated_reactions.csv` carries the kinetic estimate
`keq_kin` (with `keq_kin_basis` stating whether it is a directly reported
kinetic `K'`, a reported `Keq`, or an explicit Haldane combination of
specificity constants), the matched thermodynamic comparator `keq_thermo`, the
assay conditions, and the primary kinetic and thermodynamic sources.

## Input CSV format

One row per reversible uni--uni record. A header row is required; unknown
columns are ignored and missing optional columns may be left blank.

**Required:** `kcat_f`, `km_s`, `kcat_r`, `km_p`, and exactly one
thermodynamic comparator -- either `keq_thermo` (a dimensionless `K'`) or
`dg_thermo_kJ_per_mol` (standard transformed Gibbs energy, kJ/mol).

| Column | Meaning |
|--------|---------|
| `record_id` | identifier used in the output and summary |
| `ec_number`, `enzyme_name`, `organism`, `reaction_rhea` | provenance / normalization metadata |
| `kinetic_source`, `thermo_source` | e.g. SABIO-RK, BRENDA, TECRDB, eQuilibrator |
| `pH`, `temperature_K` | assay conditions (`temperature_K` sets `RT`; defaults to 298.15 K) |
| `kcat_f`, `kcat_r` | forward / reverse turnover (`1/s`); units cancel in the ratio |
| `km_s`, `km_p` | Michaelis constants for substrate / product (e.g. `mM`); units cancel |
| `keq_thermo` *or* `dg_thermo_kJ_per_mol` | thermodynamic comparator |
| `sigma_ln_kcat_f`, `sigma_ln_km_s`, `sigma_ln_kcat_r`, `sigma_ln_km_p`, `sigma_ln_keq_thermo` | optional 1-sigma uncertainties in natural-log units; if **all** are present, the standardized residual `z = delta / sigma_delta` is reported |

## Output columns

The output CSV repeats the identity columns and adds `keq_kin`, `keq_thermo`,
`ratio_x`, `ln_ratio`, `C_haldane`, `z`, `consistency_class`, and the
baseline measures `abs_log_error`, `squared_log_error`, `abs_ddg_kJ_per_mol`.

## Fold-discrepancy bands

From the symmetric fold error `f = max(x, 1/x)` (the manuscript fold-band table):

| Fold band | Condition | `C_Haldane` |
|-------|-----------|-------------|
| within twofold | `f <= 2` | `<= 0.25` |
| 2-5-fold | `2 < f <= 5` | `0.25 < C_Haldane <= 1.60` |
| 5-10-fold | `5 < f <= 10` | `1.60 < C_Haldane <= 4.05` |
| >10-fold | `f > 10` | `> 4.05` |

## Library usage

```python
import haldane_consistency as hc

# directly from kinetic and thermodynamic constants
keq_kin = hc.keq_kinetic(kcat_f=200, km_s=0.5, kcat_r=100, km_p=0.4)   # 1.6
score   = hc.haldane_score(keq_kin, keq_thermo=1.5)                    # ~0.002
band    = hc.classify(keq_kin / 1.5)                                    # "within twofold"

# or a full record (with optional uncertainties)
rec = hc.EnzymeRecord(record_id="R1", kcat_f=200, km_s=0.5, kcat_r=100,
                      km_p=0.4, keq_thermo=1.5)
result = hc.score_record(rec)
print(result.C_haldane, result.consistency_class)
```
