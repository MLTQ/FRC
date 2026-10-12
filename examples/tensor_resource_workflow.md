# Stage resource observations in the complete fine physical workflow

Runs the entire F64/256, B256/256 tensor baseline and verifies every saved
input/rate/face array and non-runtime physical metric against the unobserved
baseline. Adds only entry/exit wall-time and own-process peak-RSS observations
around existing current transport, current/velocity gradients, current stress,
species pressure, electromagnetic forces and metric construction.

`run` uses the existing 8 GiB / 256 MiB resource guard. Serial method bindings
restore on success or failure. A separate observation archive includes this
instrumentation and the complete physical/guard source closure. The physical
calculation retains its source identity plus this explicit observation record.

RSS is macOS cumulative process high-water memory: a rise locates an interval
that established a new maximum, but does not identify all live allocations.
Nested intervals overlap; never add their peak deltas or interpret them as
exclusive retained storage. Runtime/resource quantities are not expected to
replay bitwise. These observations target the unresolved larger-grid execution
cost; they cannot establish refinement feasibility by themselves.

No unit tests, mathematical changes, tolerance changes, numerical diffusion,
filters, atmosphere or repairs. A successful profile requires a full physical
replay; a resource interruption remains distinct from numerical rejection.
