# Lossless scientific record formatting

Writes ordinary JSON with readable object fields and compact numeric/scalar
array trees. Large limiter maps must not create a million one-number diff
lines. Lists of objects, such as case records and histories, remain indented.
The round-trip parse is checked against the original object before writing;
nonfinite values are rejected. This changes representation only, never recorded
measurements or their engine/measurement provenance. The local report records
this formatter's separate hash and formats both collected and per-case records.
