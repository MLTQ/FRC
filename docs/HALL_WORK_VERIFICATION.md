# Hall spatial verification checkpoint

**244 tests pass (84.28 s). All 42 declared spatial cases completed. The new
Hall operator is experimental and is not enabled in the compressible solver.**
Its global work identity passes, but its regular-axis consistency gate fails.
No nonlinear Hall time evolution or FRC formation is certified here.

## Independent physics references

The new finite-pressure cylindrical oracle checks all three wave branches,
both axial directions, signed bias/frequency, ideal-MHD and Cartesian parallel
limits, and approach to the incompressible limit. Physical Maxwell and original
species momentum residuals are independently differenced; true volume energy
is checked at different phases. The nonlinear polynomial/Fourier snapshot
checks variable-density Hall derivatives and the pointwise Poynting identity,
including the axis. See COMPRESSIBLE_HALL_REFERENCE and HALL_WORK_NUMERICS.

## Nonlinear refinement

All norms below retain the axis cells. Physical fields are sampled inside cells
with physical 2 pi r volume weighting; cell-mean errors are separately labeled.

| Radial cells | Current relative L2 | Physical magnetic-rate relative L2 | First-annulus g-rate relative L2 | Local energy-defect relative L2 |
|---:|---:|---:|---:|---:|
| 8 | 0.0861516 | 0.259461 | 0.202924 | 0.0261605 |
| 16 | 0.0264879 | 0.132568 | 0.211155 | 0.00689868 |
| 32 | 0.0104290 | 0.0666923 | 0.213737 | 0.00174652 |
| 64 | 0.00481440 | 0.0333967 | 0.214411 | 0.000437990 |
| 128 | 0.00235369 | 0.0167042 | 0.214581 | 0.000109583 |

Corner psi-rate maximum relative error decreases 0.0614134 → 0.0156862 →
0.00393345 → 0.000982854 → 0.000245539. The conservative g-mean volume error
decreases to 0.00586535. However, the first-annulus coefficient error tends to
**21.46%**, despite refinement and roundoff global conservation. Its volume
shrinks, allowing the global norm to conceal a persistent regular-axis error.
The experimental stencil must be revised before coupling it to fluid energy.

Local energy defect means the existing magnetic-energy rate minus the
independently evaluated conservative Hall Poynting flux rate. It would be a
numerical thermal error, not a physical heat source. It is measured, not repaired.
Its roughly second-order global decrease does not excuse the induction-axis
failure. Zero initial toroidal/poloidal perturbation and constant/stronger
carrier contrast controls are also retained; none is a replacement target.

## Work, quadrature and wave controls

Across all cases, maximum assembled Hall work divided by integral |J| |E| is
2.72060e-17. Independent archive replay of physical magnetic-energy work has
maximum backward error 2.60690e-17 against the actual absolute summation scale.
For weak waves on a strong bias, that replay scale is O(amplitude), while
Hall work is O(amplitude²); conflating them creates a misleading roundoff test.
No density floor, energy reset, Hall Joule heat or extra fluid force is present.

Reference cell quadrature orders 10/14 differ by at most 4.75452e-16 relative.
At nr=32, common-quadrature orders 3→5 change magnetic-rate error by 2.49148e-8;
orders 5→7 change it by 1.97065e-15. The axis failure is not removed by more
quadrature. Cell-constant carrier density is a separately labeled control.

For the three finite-pressure branches, positive-k physical Hall-rate errors
at nr=8/16/32 are respectively

- Branch 0: 0.159309, 0.0797043, 0.0398418.
- Branch 1: 0.140302, 0.0702920, 0.0351560.
- Branch 2: 0.160022, 0.0800576, 0.0400178.

Negative-k cases agree to numerical precision. These are instantaneous Hall
induction comparisons to finite-pressure mode fields, not propagated waves.
Their first-order physical rate accuracy remains explicit.

Direct evaluation of the quadratic Hall contribution at amplitudes
0.001/0.01/0.1 m/s gives RMS rates 3.91246115e-12, 3.91246103e-10 and
3.91246104e-8 T/s, with fitted adjacent powers 1.999999987 and 2.000000001.
The algebraically equivalent sum of opposite full-field rates loses accuracy
through cancellation: relative discrepancies are 17.93, 0.1736 and 0.001715.
Both calculations are saved and plotted; the noisy subtraction is not treated
as a physical departure from quadratic scaling.

## Reproducibility and next action

Run `examples.hall_work_measurements`, then `examples.verify_hall_work` from
the repository using the documented plotting cache environment. The frozen
production/measurement fingerprint is
`35cea1122c521983ac6808a8d8f6d27c3cf430428561e2f0de45819598a93996`.
Every configuration, case record, source identity and lossless archive hash
has been checked. Both final figures were visually inspected.

![Hall spatial error and failed axis gate](../outputs/hall_work_verification/hall_spatial_accuracy.png)

![Axis waveform, wave branches and cancellation control](../outputs/hall_work_verification/hall_axis_and_modes.png)

The next work is a regular-axis polynomial-consistent current/electric/energy
pairing, then local Hall energy transport, a stable time integrator and actual
compressible nonlinear wave evolution with time and mesh refinements.
Electron/ion thermal dynamics, finite-mass transport, physical vacuum/edge and
chamber ends, spatial coil work, kinetic ordering, converged formation and the
experimental interface remain required. The full goal stays active.
