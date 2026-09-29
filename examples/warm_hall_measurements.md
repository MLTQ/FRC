# Warm Hall protocol driver

Runs all declared pressure, boundary and coupled measurements under one frozen
transitive source identity. Optional independent groups permit concurrent
bounded runs into distinct case files; they do not publish a complete collection.
`--assemble` requires all 58 exact configurations, source identities, completion
states and NPZ hashes. Default full-protocol execution assembles only after every
case succeeds. A failed or absent case cannot be silently skipped or relabeled.

All inputs and nonlinear budgets are explicit in case JSON, units SI. Grouping
changes scheduling only. No source may change while any measurement group is
active. Report/audit/figure modules are separate from these measurement sources.
