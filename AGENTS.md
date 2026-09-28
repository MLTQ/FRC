# Agent Instructions

This project uses **bd** (beads) for issue tracking. Run `bd onboard` to get started.

## Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --status in_progress  # Claim work
bd close <id>         # Complete work
bd sync               # Sync with git
```

## Research workflow and continuation

This is a local research project with no configured remote. The initial `bd init`
template's mandatory publication workflow does not apply. Do not invent a remote
or publish research data without user authorization. Export/sync issue state
locally; record any repository-state limitation and continue scientific work.

1. Read `docs/GOAL_STATUS.md` and the relevant module companion documentation.
2. Track remaining requirements with beads; keep the overarching FRC goal open.
3. State equations, assumptions, boundaries and numerical diffusion explicitly.
4. Verify against independent analytic or manufactured solutions, not merely
   implementation-shaped assertions. Measure convergence and conservation.
5. Preserve exact vacuum outside plasma in formation runs; never quietly replace
   it with a conducting or massive numerical atmosphere.
6. Distinguish verification fixtures from physical formation experiments.
7. Update companion docs, measured results, issues and goal status after changes.
8. At a continuation boundary, leave concrete evidence and the next action.

No sub-agent delegation is required by this repository. Hardware remains abstract
prescribed waveforms. Neither closed flux nor an axisymmetric equilibrium proves
three-dimensional stability or particle confinement.
