# Supplement deposit checklist

Steps to turn this `code/` supplement into a citable, version-pinned archive and to
fill the placeholders in the manuscript's Supporting Material section
("Data and code availability").

## 1. Clean the working tree
- [ ] Remove build artifacts and scratch outputs:
      `rm -rf __pycache__ .pytest_cache scored_example.csv`
- [ ] Confirm only the intended files remain.

## 2. Put the supplement under version control
- [ ] From the repository root, `git init` (if not already a repo) and commit the `code/`
      tree together with the manuscript sources.
- [ ] Tag the release: `git tag -a v1.0 -m "Manuscript submission snapshot"`.
- [ ] Record the resolved commit: `git rev-parse --short HEAD` (short) and
      `git rev-parse HEAD` (full).

## 3. Push and archive with a DOI
- [ ] Push to a public host: `git push origin main --tags`.
- [ ] Create a Zenodo deposit: enable the GitHub--Zenodo integration and publish the `v1.0`
      release, or upload a zip of `code/` plus the manuscript. Zenodo mints a versioned DOI
      of the form `10.5281/zenodo.XXXXXXX`.

## 4. Verify reproducibility from a clean checkout
- [ ] `git clone <url> && cd <repo>/code`
- [ ] `python test_haldane_consistency.py`  ->  expect `17/17 tests passed`.
- [ ] `python score_curated.py --output /tmp/check.csv && diff curated_reactions_scored.csv /tmp/check.csv`
      ->  expect no differences (the RacE1 row scores C = 0.235).
- [ ] (Optional) `pip install -r requirements.txt && python make_real_figures.py` regenerates
      Figures 3--4 into `../figures/`.

## 5. Fill the manuscript placeholders
In the Supporting Material section, replace the placeholder tokens with the real values and
delete the `% TODO(before submission)` comment:
- [ ] `10.5281/zenodo.XXXXXXX`  ->  the minted Zenodo DOI.
- [ ] Git tag `v1.0`            ->  the actual tag, if different.
- [ ] commit `XXXXXXX`          ->  the short commit hash from step 2.
- [ ] No manuscript Lean pin to confirm: the manuscript proof is self-contained (Appendix A, Lemma 3) and does not cite or require the optional Lean material in `lean_formalization_details.md`.
