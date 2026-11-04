# hybrid3d_array_check.py

Measures the concrete float64 stencil optimization using existing physical
comparisons and actual coupled-bank initial-fill shots. No unit suite.

- `independent`: Gaussian translation with transport subcycling at32/64/128,
  and full kinetic moving-contact/constant-eta/Spitzer decay at16/32/64.
  These prescribed periodic fixtures are not physical formation experiments.
- `device`: actual coupled shot; excludes initialization and the first ten
  steps from timing. Preserves setup, endpoint snapshots, all float64 physical
  arrays, energy/loss history and source paths/hashes. Optional cProfile covers
  the same warmed interval. The default0.75us replay is not full-shot timing.

Use `.venv/bin/python -m examples.hybrid3d_array_check OUT --fused` for the
candidate. Run the same script by absolute path with `PYTHONPATH` and cwd
set to a pre-change source copy for the original control, omitting `--fused`.
The recorded closure path identifies which version actually ran. Inputs use
the absolute project output path, so this does not change the physical fill.
`--gpu 2` measures the mixed-precision backend with/without fused CPU work;
its pre-existing precision sensitivity is a separate unresolved issue.

`--dx .01 --dz .01 --steps 40 --section device` times the actual1.29million
ion startup. All synchronization and the unmodified vacuum/circuit/particle
work are included. GPU access requires authorized Metal-capable execution.
The flag is recorded and opt-in; no running source-frozen queue is changed.
Measured results and timing qualifications are in `docs/HYBRID_3D.md`.
The reference and final CPU candidate have exactly matching complete float64
physical arrays at both coarse and fine endpoints. A repeated fine reference
exposes timing variability; it is retained rather than selecting one result.
