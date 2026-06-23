# Lean formalization details

This note identifies optional, public Lean 4 material that independently checks the
d'Alembert/cosh classification route used by the reciprocal-cost uniqueness theorem
(Theorem 2). The manuscript proof is self-contained --- the classification is stated and
proved as Lemma 3 in Appendix A --- so this formalization is neither cited by nor required
by the manuscript; it is provided here only as a supplementary independent check.

## Repository and pinned state

- Repository: https://github.com/jonwashburn/shape-of-logic
- Pinned commit: `edcd65cc4b88db0cb23c5918bfcfbee45fcbbc1a`
- Short commit: `edcd65c`

## Exact checkout and build commands

```
git checkout edcd65cc4b88db0cb23c5918bfcfbee45fcbbc1a
lake build IndisputableMonolith.Foundation.DAlembert.Inevitability
```

## Relevant module and key declarations

- Module: `IndisputableMonolith.Foundation.DAlembert.Inevitability`
- Source file: `IndisputableMonolith/Foundation/DAlembert/Inevitability.lean`

Representative declarations in the pinned commit:

- `bilinear_family_forced`
- `bilinear_family_reduction`
- `axiom_bundle_necessary`

## Related uniqueness module (supporting context)

The reciprocal-cost uniqueness statement itself is also recorded in a companion
module at the same commit:

- Module: `IndisputableMonolith.CostUniqueness`
- Source file: `IndisputableMonolith/CostUniqueness.lean`
- Representative declarations: `T5_uniqueness_complete`,
  `unique_cost_on_pos_from_rcl`, `unique_cost_on_pos`, `Jcost_regularity_cert`

Both modules build at the pinned commit. The manuscript does not cite or rely on this
material; the `Foundation.DAlembert.Inevitability` module corresponds to the d'Alembert/cosh
classification route, with `CostUniqueness` as the companion uniqueness statement.
