# Clear memory-legend placement

Preserves the original scientific report and its artifacts. A separately
archived renderer moves only the memory panel's legend above the axes, clear
of the8GiB guard line and numeric bar labels. Original physical inputs, metrics,
captions, curves, bars and all acceptance requirements stay unchanged.

`run` temporarily scopes the legend placement and restores bindings on exit.
Its source archive includes both renderers; visually inspect the new PNG.
This layout change does not introduce scientific verification or unit tests.
