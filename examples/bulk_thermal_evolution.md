# Compression and full-transit transport protocol

Ten coupled density/two-entropy runs cover compression and full axial translation
at nr=8/16/32, plus half time step, quadrature, reverse compression and degree 3.
The prescribed velocity acts on the current numerical density, never an analytic
replacement. All accepted densities and species pressures validate as positive.

Each accepted step retains particles, physical entropy, internal energy, pressure
extrema and the independently accumulated pressure-work ledger. Five saved
states receive independent physical continuum norms at two quadratures. Initial
and evolved errors remain distinct. Last pre-step state permits an exact replay.
The velocity is prescribed; no momentum/total-plasma energy conservation or
magnetic formation is implied. Failures and their tracebacks remain recorded.
The corrected driver also saves the last accepted state/history, attempted step,
original exception and partial snapshots if a numerical step fails. Earlier
baseline failures predate this retention and contain their original traceback;
that limitation is explicit, not backfilled with invented accepted arrays.
The command runs up to two independent cases in separate processes (selectable
with --workers), saves each case before collection, and emits the final protocol
in its fixed index order. No state is shared between numerical trajectories;
this changes orchestration only. All jobs must finish before the source guard
and complete audit can pass.
