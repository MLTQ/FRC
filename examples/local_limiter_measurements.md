# Local-limiter comparison measurements

Runs the unchanged independent cylindrical shock through nr=384, same-source
global controls at 96/192, half-CFL at 192, doubled gas half-length at 96, and
the disabled-limiter failure. Identical driven pulses and independent smooth
solutions check that inactive limiting retains their accuracy and source budgets.
Uses the common complete input/history/error and sparse limiter-map recorder.

Each case is immediately persisted in its own file. Reuse requires the package,
common measurement helper, this case-list source and environment fingerprints
to match. `--case` selects named cases so independent expensive runs can use
separate processes without writing the same file. Historical conservative and
split-study files are never overwritten. Failed runs remain accepted-time
records and cannot count as final-time convergence evidence.
The recorded nr=384 case required about 36 minutes including setup on this
CPU environment. Use the ordinary test suite for a fast check and `--reuse`
for unchanged measurements; do not mistake long execution for completion.
