# Release Verification for v1.3

This file records the local checks performed before publishing the synchronized
GitHub/Zenodo v1.3 release for:

**A Reciprocal Reporting Scale for Haldane Consistency in Reversible Enzyme
Kinetics: A Curated Proof of Concept**

## Completed Local Checks

- `curated_reactions.csv` excludes the provenance-untraced `PGI_88TEW` record.
- `curated_reactions_scored.csv` was regenerated from `curated_reactions.csv`.
- Fold-discrepancy labels use the manuscript bands: `within twofold`,
  `2-5-fold`, `5-10-fold`, and `>10-fold`.
- Generated code outputs no longer report the legacy standardized-score field;
  the optional uncertainty output is the standardized residual `z`.
- `README.md` and `CITATION.cff` use the manuscript title.
- `CITATION.cff` is set to `version: "1.3"`.
- Supplementary Tables S1-S3 are included as explicit CSV files.
- Expected SHA-256 checksums are recorded in `EXPECTED_OUTPUT_HASHES.sha256`.

## Commands Run

```bash
python score_curated.py --output curated_reactions_scored.csv > score_curated_output.txt
python make_real_figures.py > make_real_figures_output.txt
python test_haldane_consistency.py
```

The test run reported `17/17 tests passed`.

## Checksum Verification

From the `code/` directory:

```bash
shasum -a 256 -c EXPECTED_OUTPUT_HASHES.sha256
```

Expected result: every listed output reports `OK`.

## Post-Upload Checks

- Create GitHub tag/release `v1.3` from the uploaded synchronized files.
- Archive the same `v1.3` release on Zenodo.
- Replace manuscript DOI placeholders with the minted Zenodo version DOI.
- Confirm the concept DOI still resolves to the latest release.
