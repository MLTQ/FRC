# Explicit historical global-limiter driver

The original conservative report describes global convex limiting. After the
solver default changes to local fallback, rerunning that report must still
select the global method. This driver reads the tracked historical study's
complete input records as a manifest and explicitly sets `flux_limiter=global`
for every case, including inactive/disabled controls.

The manifest is loaded before any requested output is written, avoiding a
partially overwritten protocol. Only configuration values reach the existing
validated measurement factory; no code is loaded from the records. Driver,
normalized protocol, common measurement and engine fingerprints are recorded.
Reuse requires an exact match. A named subset can be measured for a smoke check.
The default local study has its own driver and never invokes this one.
