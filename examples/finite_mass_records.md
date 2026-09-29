# Conservative finite-mass provenance

Recursively hashes transitive local imports from the coupled integrator and
rate and trajectory measurement entry points, plus this provenance implementation. Records
Python/NumPy/SciPy versions. It is separate from all previous frozen current,
thermal, Hall and moving-rate studies; no prior source identity is redefined.
Measurements assert unchanged source at start/end. Report-only helper/test
sources are included explicitly when used. The initial rate-only snapshot has
its own earlier identity; it is never relabeled when driver coverage expands.
