# Current simulation field display — 2026-10-12

Open `outputs/trajectory_display/index.html` in a WebGL2 browser for the offline
interactive display. The standalone scientific image is
`outputs/trajectory_display/current_state_3d.png` (2880 × 1872), visually inspected.
The new files leave both earlier showcases and the user's poster unchanged.

The viewer shows four actual saved states: 0, 16.6667, 25 and 28.3333 µs from
the refined F8/64, B32/64 contact trajectory. The last is step 17; step 18
rejected and the requested 100 µs duration was not reached. This is a filled
periodic cylinder, radius 35 mm, length 240 mm, with a 5 mT axial bias.
These verification dimensions are not a proposed laboratory specification.

The left interactive panel revolves the magnitude of all three weak-Ampere
current components. The right follows the full B vector, including Btheta.
Common linear color scales cover all times. The optional per-frame current
scale exposes very weak early currents and explicitly updates its legend.
Camera rotation, cutaway, opacity and PNG export affect presentation only.
The static PNG uses colored cutaway surfaces instead of optical volume rendering.

## Measured display samples

| Time (µs) | Sampled peak current (A/m²) | Sampled total B range (mT) |
|---:|---:|---:|
| 0 | 1.85680 | 4.9999904–5.0000237 |
| 16.6667 | 18.7345 | 4.9999879–5.0000273 |
| 25 | 42479.3 | 4.9969708–5.0029488 |
| 28.3333 | 835178 | 4.9288188–5.0711423 |

Final sampled minimum Bz is +4.9288187 mT. The bias dominates the magnetic
geometry, so the displayed lines are nearly axial. There is no demonstrated
field reversal or FRC formation. Large current and a replayed failure are not
evidence of physical instability, experimental suitability or 3D stability.

## Verification and presentation limits

The original full trajectory and rejection replay are hash-checked before
sampling. Each current uses the unchanged production weak Ampere operator.
Independent degree-five root-polynomial reconstruction agrees with production
quadrature to relative vector norms ≤5.960e-17 for B and ≤5.937e-15 for J.
All four float64 display NPZ files, complete metrics and compact JSON replay
byte identically in a fresh process. See `trajectory_display_data`,
`trajectory_display_data_replay` and its `_verification` output.

The 256 × 512 display grid and linear interpolation are approximations. On
the declared interleaved off-grid points, final vector relative differences
are 0.00117959 for B and 0.0330405 for J. Those are unweighted sample norms,
not cylindrical physical-error norms. Sampled peaks/ranges are not exact
extrema; production quadrature samples different positions. Thin outer
half-pixel strips clamp to the nearest sample in the interactive texture.
The 400/800-step RK4 line-position difference is ≤3.007e-10 m over this
interpolated field; that measures tracing error only. No physical state was
filtered, replaced or advanced.

The PNG was rendered end-to-end and inspected. Interactive browser inspection
could not run while the Mac was locked; no UI verification is claimed.
The local viewer is served on 127.0.0.1:8786, session 74537, for the requested
display. The HTML also works offline; only direct saving beside the HTML
requires the optional local server. No unit tests were added or run.

The first data computation completed and wrote a terminal summary; its outer
`time -l` wrapper returned 1 because macOS denied a clock query. The fresh
Python-only complete replay and replay verifier both exited 0. That wrapper
failure is not a failed scientific reconstruction or an observed resource peak.
