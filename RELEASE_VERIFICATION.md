# Release Verification for v1.4

This file records the local checks performed before publishing the synchronized
GitHub/Zenodo v1.4 release for:

**Haldane Consistency as a Reciprocal Reporting Scale for Reversible Enzyme
Kinetics: A Curated Benchmark**

## Completed Local Checks

- `curated_reactions.csv` excludes the provenance-untraced `PGI_88TEW` record.
- `curated_reactions_scored.csv` was regenerated from `curated_reactions.csv`.
- Fold-discrepancy labels use the manuscript bands: `within twofold`,
  `2-5-fold`, `5-10-fold`, and `>10-fold`.
- Generated code outputs report the standardized residual `z` for the optional
  uncertainty output (no legacy standardized-score field).
- `README.md` and `CITATION.cff` use the combined-manuscript title.
- `CITATION.cff` is set to `version: "1.4"`.
- Supplementary Tables S1-S3 are included as explicit CSV files.
- Expected SHA-256 checksums are recorded in `EXPECTED_OUTPUT_HASHES.sha256`.

## Commands Run

```bash
python score_curated.py --output curated_reactions_scored.csv > score_curated_output.txt
SOURCE_DATE_EPOCH=1700000000 python make_real_figures.py > make_real_figures_output.txt
python equilibrator_estimates.py > equilibrator_estimates_output.txt   # optional; needs equilibrator-api
python test_haldane_consistency.py
```

The test run reported `20/20 tests passed`.

## Checksum Verification

From the `code/` directory:

```bash
shasum -a 256 -c EXPECTED_OUTPUT_HASHES.sha256
```

Expected result: every listed output reports `OK`. The two figure PDFs are
byte-reproducible only when `SOURCE_DATE_EPOCH` is pinned as shown above.

## Post-Upload Checks

- GitHub tag/release `v1.4` published:
  https://github.com/megansimons/haldane-consistency-supplement/tree/v1.4
- Archived on Zenodo under version DOI 10.5281/zenodo.21084024; the concept DOI
  10.5281/zenodo.20790110 resolves to the latest release.
- The v1.4 version DOI is recorded in `CITATION.cff` and the manuscript.
