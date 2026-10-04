# Species convergence report

Reads completed immutable cases and their terminal audit, checks raw hashes,
and plots six scientific panels without changing data. Reports original rate,
pressure, axis/wall/max and full-transit signal error, preserves the inaccurate
point-transfer control, and separates grid from continuous quadratic drift.
Independent quadrature comparison excludes sampled maxima because the sampling
locations change between orders. Three successive endpoint differences measure
time order on the same spatial mesh. The report source hash is separate from
the frozen calculation source. Existing reports cannot be overwritten.
Constant-state backward error is recomputed from the saved weighted matrix as
norm(A*sqrt(V*n))/(norm(A)*norm(sqrt(V*n))). The control record's older
`constant_rate_scaled` omits the vector norm and carries sqrt(particle) units;
it is retained as raw data and is not presented as a relative backward error.
