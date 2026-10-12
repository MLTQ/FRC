# Release expired current-transport workspace

`ReleasedCurrentTransport` inherits original geometry/differential/face setup.
Its pair and advection methods reproduce every original arithmetic expression
and call order, adding only del statements after the last use of gradient,
stress, sampled current, velocity, drift and face temporaries. Returned arrays
and all physical terms retain their original meanings and definitions.

The source AST must agree with CurrentTransport after removing these del
statements. Full signed original/new physical workflows must then agree
bitwise and fresh candidate replays must pass. The source check alone does
not establish numerical or physical accuracy. Measure actual peak memory;
allocator behavior can retain released storage.

No external state or argument is mutated by the added deletions. All modes,
physical masses/density, central face terms, cylindrical connection terms,
work adjoints, residual checks and budgets remain unchanged. No diffusion,
filter, artificial atmosphere or repair. Original source remains frozen.
This is opt-in workspace management, not a new transport model or time method.
