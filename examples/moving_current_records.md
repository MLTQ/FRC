# Frozen moving-current provenance

Traverses actual local imports from the measurement/reference/production entry
points, hashes every transitive Python source, and records Python/NumPy/SciPy
versions. Kept separate from the historical current_thermal records so their
frozen identity and replay remain unchanged. Measurements require identical
provenance at start and finish; report-only helpers carry separate hashes.
The finite_mass_moments continuum contract is an explicit additional root.
