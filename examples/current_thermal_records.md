# Current/thermal records and source identity

Traverses transitive local import/from dependencies, including function-local
imports, from the actual driver. Every source SHA256 and runtime version are
recorded. Sources must be unchanged at every save. Full configurations,
unrounded accepted histories, explicit failure status and compressed lossless
arrays carry their own hashes. Rejected trials are never labeled completed.
Report and test sources are hashed separately; they do not alter measurements.
