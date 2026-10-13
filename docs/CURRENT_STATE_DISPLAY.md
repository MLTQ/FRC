# Updated 3D field display — 2026-10-12

The current entry points are [the offline viewer](../outputs/current_state_display_v2/index.html)
and [the inspected scientific PNG](../outputs/current_state_display_v2/current_state_3d.png).
Both show the last accepted trajectory state: step 17, 28.3333 µs. The newer
spatial-contribution and local-current audits investigate that trajectory;
they have not advanced it to a later time.

The viewer retains four real saved states, shared physical color scales, full
magnetic field lines including the 5 mT bias, and all-component current magnitude.
New side/end/perspective presets and keyboard rotation/zoom supplement the
existing linked cameras, cutaway, opacity, per-frame scaling and PNG export.
The latest audit status is now on the page and standalone PNG. The offline HTML
contains all data and scripts; it makes no external data requests.

At the final displayed time, sampled peak |J| is 835.178 kA/m² and sampled
minimum Bz is +4.92882 mT. The full B magnitude spans 4.92882–5.07114 mT.
These are display-sampled extrema. Magnetic deflection is not exaggerated.
The filled periodic verification cylinder is 70 mm in diameter and 240 mm long;
these are verification dimensions, not a laboratory design recommendation.

Step 18 failed to converge, so the requested 100 µs trajectory remains incomplete.
The axial field stays positive and no FRC formation is demonstrated. The latest
fully replayed audits expose internal current divergence and radial normal-current
jumps in the saved error direction; their role in rapid growth remains unresolved.
See [the scientific audit](SPATIAL_PARTS_VERIFICATION.md). Normalized diagnostic
directions are not substituted for evolved fields in the display.

## Verification

Both complete display builds ran in fresh processes. HTML, PNG and manifests
reproduce byte-identically. Every direct manifest input and the preceding
checkpoint's scientific outputs were hash-checked. The final 2880 × 1872 PNG
was visually inspected. The first closer-camera PNG clipped cylinder ends;
that source/output is preserved separately, and the delivered v2 corrects only
the view zoom. The generated interactive HTML is identical between the two
presentation variants.

Assembled JavaScript passes Node syntax checking. Browser interaction could not
be checked because the Mac was locked; syntax checking and static rendering
do not verify the interactive UI. The PNG preview was queued in Codex. No unit
tests were added or run. No solver or accepted state was changed.

The independently verified sampling, interpolation errors, cell-center edge
clamping and line tracing limits remain exactly those in
[TRAJECTORY_DISPLAY.md](TRAJECTORY_DISPLAY.md). The current view uses approximate
cutaway surfaces in the PNG and an illustrative optical volume in WebGL.
Neither is a measurement of optical emission or a 3D physical simulation.

Rebuild to a new destination with:

```sh
PYTHONDONTWRITEBYTECODE=1 VECLIB_MAXIMUM_THREADS=1 \
MPLCONFIGDIR=/private/tmp/plasma-bench-mpl MPLBACKEND=Agg \
.venv/bin/python -m examples.current_state_display_v2 --output outputs/my_field_display
```

Completed output directories are protected from overwrite. Add `--serve` for
an optional localhost server on port 8787, only if that port is free. The HTML
works offline without the server. The full FRC goal remains open; the next
scientific action is the independent weak species-entropy interface budget in
[SPECIES_INTERFACE_DIAGNOSTIC_DESIGN.md](SPECIES_INTERFACE_DIAGNOSTIC_DESIGN.md).
