# Original user project request

Preserved verbatim from the original attached request for use by the next thread.
The later active objective and verified project status are recorded in
GOAL_STATUS.md.

---

We are starting a new software/research project: an interactive “virtual plasma bench” for exploratory plasma-formation and magnetic-topology experiments.

BACKGROUND / MOTIVATION

I have a home plasma/vacuum lab and have been thinking about minimal-lab approaches to forming FRC-like (field-reversed configuration) plasmas. One conceptual experiment is:

- Long cylindrical/quartz vacuum chamber.
- Low-pressure gas fill (argon initially, but species should ultimately be selectable).
- Pre-ionized plasma.
- Approximately uniform axial bias field, conceptually from a Helmholtz/solenoidal magnet arrangement.
- One or more axisymmetric pulsed coils around the chamber.
- Fast changing external fields induce toroidal plasma current J_theta.
- We want to determine whether that induced current can reverse the axial magnetic field and produce closed poloidal flux / an FRC-like configuration.

We also discussed comparing conventional theta-pinch-like unipolar pulses against oscillatory bursts resembling the magnetic waveform produced by a resonant ZVS induction-heater circuit. For now, DO NOT design hazardous pulsed-power hardware. Treat external hardware as prescribed current or magnetic-field waveforms I(t), and determine what waveform requirements the plasma physics implies.

The larger goal is NOT merely an FRC simulator. I want a general-purpose interactive plasma sandbox where I can define a chamber, gas/plasma, coils, fields, and waveforms, hit Run, and watch the resulting electromagnetic/plasma dynamics.

Eventually I want to experiment with theta pinches, FRCs, mirrors, cusps, spheromak-ish configurations, multipoles, RF fields, strange coil arrangements, and eventually fully 3-D configurations.

CORE DESIGN PRINCIPLE

Use progressive physics.

Do not attempt a gigantic kinetic simulation immediately. Build increasingly sophisticated models behind clean interfaces:

Stage 0 — Axisymmetric vacuum electromagnetics.
Stage 1 — Conducting/resistive plasma + magnetic diffusion / induced J_theta.
Stage 2 — Dynamic resistive MHD.
Stage 3 — Energy equation + realistic transport + ionization/neutrals/radiation.
Stage 4 — Hall/two-fluid effects where appropriate.
Stage 5 — eventually 3-D and/or kinetic/hybrid extensions.

The program should always tell us when the assumptions of the currently selected model are questionable.

FIRST TECHNICAL SPECIFICATION

Use an axisymmetric cylindrical (r,z) domain.

Represent the poloidal magnetic field using poloidal flux psi:

    B_r = -(1/r) dpsi/dz
    B_z =  (1/r) dpsi/dr

(up to whichever precise psi convention we adopt — document it carefully).

This representation should preserve div B = 0 by construction.

Flux contours should directly visualize magnetic surfaces, O-points, X-points, separatrices, and closed-flux regions.

Separate the architecture approximately as:

    geometry
       ↓
    coils + waveforms
       ↓
    vacuum electromagnetic field
       ↓
    plasma model / solver
       ↓
    diagnostics
       ↓
    visualization / parameter sweeps

The coil/electromagnetic subsystem should NOT be tightly coupled to a particular plasma solver.

STAGE 0: VACUUM ELECTROMAGNETICS

Implement arbitrary circular axisymmetric coils.

Each coil should have parameters such as:

- radius
- axial position
- turns / effective ampere-turns
- polarity
- current waveform I(t)

Support waveform objects such as:

- DC
- Gaussian-ish pulse
- finite-rise unipolar pulse
- RLC/damped-sinusoid pulse
- sinusoidal burst
- arbitrary sampled waveform
- arbitrary user callback/function if practical

Calculate:

- psi(r,z,t)
- B_r
- B_z
- magnetic-field magnitude
- dpsi/dt
- induced azimuthal electric field E_theta

Use analytic circular-loop expressions involving elliptic integrals where appropriate rather than approximating every coil as a pile of discrete wire elements.

Handle r=0 correctly.

Add unit tests against known limiting cases, especially the on-axis field of a circular current loop.

STAGE 1: CONDUCTING PLASMA

Add a simple prescribed plasma region with finite conductivity.

Initially it is acceptable to use a configurable conductivity sigma, before introducing detailed transport models.

Evolve induced toroidal current and magnetic diffusion consistently.

The central experiment is:

Can an external pulse induce enough J_theta to make B_z reverse locally?

Track:

- J_theta(r,z,t)
- total B_z
- minimum on-axis B_z
- B_z / B_bias
- magnetic energy
- Ohmic heating
- diffusion time
- formation/disappearance of closed psi contours

We should be able to compare:

A) fast unipolar theta-pinch-like pulse
B) underdamped RLC pulse
C) several-cycle sinusoidal/ZVS-like burst

Prefer comparisons normalized in physically meaningful ways, e.g. equal peak external field, equal delivered magnetic energy where applicable, etc.

STAGE 2: RESISTIVE MHD

Then implement axisymmetric compressible resistive MHD sufficiently well to allow:

- J x B acceleration
- radial compression
- axial motion
- pressure evolution
- magnetic advection
- resistive magnetic diffusion

Do NOT silently add artificial physics simply to make simulations look stable.

Clearly document numerical diffusion and stabilization.

Strongly favor conservative formulations and verification problems.

If implementing the full solver immediately is too large, scaffold the interface and build it incrementally.

PHYSICS / VALIDITY DIAGNOSTICS

Every simulation should calculate and display relevant characteristic scales and dimensionless parameters whenever possible:

- plasma beta
- Alfven speed
- Alfven transit time
- magnetic Reynolds number
- Lundquist number
- resistive diffusion time
- ion/electron gyrofrequencies and gyroradii where parameters permit
- collision frequencies / mean free paths once appropriate
- Hall parameter / ion skin depth when we approach two-fluid regimes
- grid spacing relative to important physical scales

The program should WARN when the selected model is likely invalid.

For example, if ion skin depth becomes comparable to the simulated structure size, tell us that single-fluid resistive MHD may be questionable.

DIAGNOSTICS / TOPOLOGY

Build topology diagnostics as first-class components.

Automatically detect where possible:

- B_z reversal
- magnetic O-points
- X-points
- separatrices
- closed psi regions
- approximate separatrix radius/volume
- lifetime of closed flux
- magnetic-axis location
- peak J_theta
- excluded-flux radius
- magnetic energy
- plasma kinetic/internal energy when those models exist
- energy conservation/error budget

Do not equate “field reversal” with successful FRC formation. Report the individual criteria separately.

VISUALIZATION

I eventually want this to feel like a virtual laboratory.

For the first implementation, visualization can be Python/Matplotlib or another pragmatic stack.

Useful views:

1. r-z cross section colored by B_z with psi contours.
2. J_theta heatmap.
3. density.
4. pressure/temperature.
5. velocity vectors once dynamic MHD exists.
6. time traces for B_z at virtual probes.
7. energy histories.
8. separatrix/O-point/X-point tracking.

Eventually I want an interactive UI where coils can literally be moved around a virtual chamber and their waveforms edited, but DO NOT let UI development delay the physics engine.

Everything should remain scriptable.

PARAMETER SWEEPS

Design configuration objects so we can automatically sweep:

- bias B
- plasma density
- gas pressure
- initial ionization
- conductivity/resistivity
- pulse amplitude
- pulse rise time
- pulse duration
- oscillation frequency
- damping
- coil radius
- coil axial position
- spacing between formation coils
- chamber dimensions

Produce phase diagrams showing regions where:

- no reversal occurs
- transient B reversal occurs
- closed flux occurs
- closed flux persists for some specified duration

SOFTWARE ENGINEERING

Please create a clean repository/project rather than one giant notebook.

Suggested conceptual layout:

    plasma_bench/
        geometry/
        coils/
        waveforms/
        electromagnetics/
        plasma/
        diagnostics/
        numerics/
        visualization/
        configs/
        tests/
        examples/

But change this if you have a better architecture.

Use SI units internally.

Use dataclasses/Pydantic/etc. for explicit configuration rather than mystery dictionaries.

Make simulations reproducible.

Include metadata with output files describing the complete configuration.

Add tests continuously.

Where an analytic solution exists, use it as a regression test.

Performance matters because this is intended as an exploratory sandbox. Vectorize/JIT where useful, but prioritize correctness first.

JAX, NumPy/SciPy, Numba, PETSc, finite-volume libraries, etc. are all options. Investigate what makes sense before committing.

GPU acceleration is attractive eventually; the development machine has an RTX 4090.

IMPORTANT RESEARCH ATTITUDE

This is an exploratory scientific codebase.

Do not make a pretty animation and declare success.

At each stage:

1. State the governing equations.
2. State assumptions.
3. State boundary conditions.
4. State discretization.
5. Derive or cite expected behavior.
6. Implement verification tests.
7. Quantify numerical error.
8. Only then run the interesting plasma experiment.

Distinguish carefully between:

- magnetic-field reversal
- closed magnetic flux
- an MHD equilibrium
- an FRC-like transient
- a stable FRC
- particle confinement

Those are NOT interchangeable.

FIRST TASKS

Start by inspecting the environment and creating the project.

Then:

1. Write docs/PHYSICS.md defining the staged model and conventions.
2. Write docs/ARCHITECTURE.md.
3. Decide on the numerical/software stack and justify the decision.
4. Implement Stage 0 completely.
5. Build circular-loop field calculations with verification tests.
6. Implement waveform objects.
7. Create a simple example consisting of:
       cylindrical chamber
       uniform axial bias field
       two formation coils
8. Visualize psi and B_z through a pulse.
9. Verify the vacuum solution numerically.
10. Begin Stage 1 only after Stage 0 tests pass.
11. Implement the simplest defensible conducting-plasma magnetic-diffusion model.
12. Run a small parameter comparison between:
       unipolar pulse
       damped RLC pulse
       sinusoidal/ZVS-like burst
13. Report what the model actually says rather than trying to force FRC formation.

Create a TODO / research-questions document for unresolved physics questions instead of making arbitrary choices.

Some obvious open questions include:

- What plasma conductivity model is appropriate during formation?
- How should pre-ionization be represented?
- At what point is Hall MHD required?
- How important are end effects?
- Which boundary conditions best represent a quartz tube surrounded by vacuum?
- When does a fluid treatment cease to be credible at the low densities/pressures we want to explore?
- What objective definition of “FRC formation” should the automated detector use?
- How should equal-energy comparisons between radically different coil waveforms be defined?
- What quantities would actual B-dot probes, Rogowski coils, photodiodes, spectroscopy, etc. measure so the simulator can eventually compare directly with experiment?

Keep those uncertainties explicit.

Finally, give me a concise summary of what was implemented, what verification tests pass, what assumptions remain, and what the next physics milestone should be.

Do not design or provide construction instructions for a hazardous high-energy pulsed-power supply. We are currently determining the desired coil-current / magnetic-field waveforms in simulation; hardware engineering can remain abstract.
