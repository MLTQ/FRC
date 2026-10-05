# Complete diagnostic replay verification

`replay(kind, original, repeated, destination)` executes the complete wall or
independent-continuum diagnostic into a fresh directory, then compares every raw
array and non-runtime metric. `verify` can certify already separately repeated
outputs without rerunning expensive identical work; the CLI uses this form.
Both complete summaries, source/runtime identities, current original source
hashes and all raw-file sets must match. Array shape/value equality is checked
for every key, with original and repeated raw hashes saved.

Writes a separate immutable verification source archive and concise success
summary only after all checks pass. Existing calculations and copies remain
untouched. Failure never writes a successful verification record. Comparison
certifies numerical replay, not physical accuracy or laboratory suitability.
