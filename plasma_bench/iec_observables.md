# iec_observables.py

Full-population read-only reconstruction for saved IEC/HOPS frames. No solver
rerun, charge modification, field smoothing or new symmetry constraint.

- `ViewMesh` builds geometry, native P1 gradients, face adjacency and lumped
  node volumes Vi=sum(tet volume/4). It does not factor or evolve Poisson.
- `moments` deposits ns=sum(wp Ni)/Vi and Js=sum(qs wp vs Ni)/Vi from all saved
  particles of each species. Volume-integrated number and current equal the
  independent particle sums to roundoff. This velocity moment is not a
  charge-conserving trajectory deposition and does not claim discrete continuity.
- `snapshot` differentiates saved phi in the native FEM space, traces E lines,
  and samples current/density/charge on labelled display grids. Current glyphs
  have normalized display lengths; color carries magnitude. Slice y=0 is 41²;
  neither grid establishes physical resolution or particle convergence.
- Phase space is (distance from measured ion centroid, relative radial speed)
  for every live ion. Absorption biases the population; these are not tracks.
  Radial density profiles are explicitly spherical shell averages about the
  chamber origin, not an assumption used in evolution.
- Optional magnetic output is the regularized plasma-only Biot–Savart
  contribution in `iec_view_fields`. No total/circuit/steel field or EM force.
- `ObservableStore` accepts local run identifiers and valid saved-state indices
  only. It bounds two mesh/eight frame caches and serializes expensive work.
  It never writes scientific states or launches/stops jobs. The HTTP adapter
  retains loopback host checks; browser requests read only local saved data.
