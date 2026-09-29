# Warm Hall source and field archives

Provenance follows transitive local `from` imports from the measurement driver,
including imports inside functions, and records every resulting source SHA256
plus Python/NumPy/SciPy versions. Full case configurations and lossless NPZ
field hashes are strict JSON. Sources are checked before each save. Report and
plotting sources are separately hashed after the engine protocol is frozen.

Every accepted step is retained in compact history columns. No data are rounded
for serialization and no failed case is marked complete. A completed record
means the declared measurement finished, not the full FRC goal is achieved.
