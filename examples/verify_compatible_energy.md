# Compatible energy report driver

Requires all 25 configured cases and four completed constant-bias time audits:
nr=24 at CFL .2/.1/.05 and nr=192 at .2. The coarse actual MHD case tests
temporal refinement separately from the strongly limited fine boundary layer.
Uses exact source/driver/environment fingerprints for reuse.
Two separately fingerprinted first-order spatial controls at nr=192 test the
same initial data at CFL .2/.1 and expose the baseline's spatial accuracy cost.
The historical local study is retained with its own provenance. The report
separates unchanged-shock profile convergence, absolute and transfer-relative
energy error, limiter/time/boundary controls, independently predicted RK defects,
identical source pulses and analytic smooth solutions. Failed cases remain
visible and are not interpreted as final-time convergence points.

Writes complete compact metrics, a standalone scientific plot and a Markdown
report. Renderer and serialization hashes are separate from numerical hashes.
The plot module owns layout and must receive visual QA. Source-normalized
energy error cannot establish resolved driven flow; semidiscrete compatibility
cannot establish temporal convergence, moving interfaces or formation validity.
Locality tables use actual accepted Euler-stage counts, maximum affected cell
fractions, fallback passes and direct replacements in the independent error window.
Six separately recorded boundary-strip probes retain the second-order interior
and report forced order reduction independently of additional adaptive limiting.
They are explicitly experimental pending driven-source and finer-boundary tests.
Reuse requires a final diagnostics record and matching engine/probe fingerprints.
Reported minimum densities retain their SI units and are never floors.
