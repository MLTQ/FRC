# Latest field poster with unclipped geometry

## Purpose
Presentation-only successor to `current_state_poster.py`. The first rendering's
1.17 box-aspect zoom clipped cylinder ends at the axes boundary. This variant
uses the original zoom of 1.0. The original source and clipped output remain
preserved; no physical input or numerical reconstruction changes.

## Components and contracts
`render(source, summary, data, target)` retains the complete contract documented
in `current_state_poster.md`: input verification belongs to the assembler;
current uses all three components, magnetic paths use full B and actual scale,
and the saved time is 28.3333 microseconds. Common linear scales, dimensions,
step-18 rejection, audit status and verification scope remain visible.

The only plotting change is the view zoom. End-to-end generation, PNG inspection
and an independent process replay are required before reporting completion.
No unit tests. Preserve this source after executing it.
