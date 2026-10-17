# Python package configuration

`pyproject.toml` defines the local plasma-bench package, its command-line entry
point, Python requirement and numerical/plotting dependencies. Research drivers
also run as modules from the repository. No external service or publication is
required for a simulation.

The optional `precision` extra pins the exercised gmpy2 2.3.2 backend. Install
with `python -m pip install '.[precision]'`, then select it before process start
using `FRC_GEOMETRY_BACKEND=mpfr`. Standard-library Decimal remains the default;
precision selection and dependency versions are recorded in each calculation's
provenance. Both use float64 physical coefficients and precise coordinates.

The existing `dev` extra and pytest configuration are historical tooling. The
current user instruction is to verify through complete simulations and
independent comparisons, without adding or running unit tests.
