"""Assemble completed resistive profiles, moving cases, controls and figures."""
from pathlib import Path
import hashlib
import matplotlib
from .resistive_analysis import read_study,comparisons
from .resistive_plot import error_figure,field_figure
from .verification_records import write_record


def report(decay,rows,controls):
    lines=["# Resistive moving-plasma verification","",
           "The moving material/vacuum reference now includes scalar physical resistivity, consistent weak current, positive per-cell Ohmic heating and atomic Strang composition with ideal motion. Adaptive midpoint substeps control magnetic-field, current and local-heat error. No exterior fluid, pressure/density floor, numerical diffusivity or energy repair is introduced. The [numerical derivation](RESISTIVE_MOVING_NUMERICS.md) gives the equations, source-work contract and limits.","",
           "All 13 stationary diffusion cases and 20 moving cases are complete and source/configuration/field-hash checked. These are unit-scale mathematical verification inputs with rho=1 kg/m³ and gamma=5/3, without a chosen species/temperature/transport regime. The Stage 2 formation entry remains gated. Older ideal and fixed-chamber studies retain historical sources.","",
           "## Independent continuum diffusion, current and heat","",
           "A fixed conducting sphere of R=1 m has D=eta/mu0=0.1 m²/s and initial central a=A=sqrt(mu0) T. For k=n pi/R, a=3A j1(ks)/(ks) exp(-Dk²t) inside, with a continuously matched vacuum dipole outside. Jtheta=k² r a/mu0 inside and zero outside. Full-space energy is 6pi A²R/(mu0 k²) exp(-2Dk²t); integrated heat is its analytic loss. Independent Maxwell/Faraday/Ohm, interface and all-space quadrature tests validate this reference, separately from the numerical operator. See the [reference derivation](MOVING_RESISTIVE_PLAN.md).","",
           "The geometry is explicitly held fixed for this component test. Its magnetic force is not balanced by fluid pressure; this is not a freely evolving equilibrium. B is evaluated from actual Q1 gradients, recovered current from the consistent weighted mass solve, and local heat from integrated eta J². At t=0.1 s for the first radial mode:","",
           "| Radial cells | B volume L2 | Flux volume L2 | Current volume L2 | Local heat volume L2 | Boundary-current scaled max error | Axis-field relative error |",
           "|---:|---:|---:|---:|---:|---:|---:|"]
    mesh=[r for r in decay["cases"] if r["name"].startswith("mesh_")]
    keys=("field_relative_volume_L2","flux_relative_volume_L2","current_relative_volume_L2","local_heat_relative_volume_L2",
          "maximum_boundary_current_scaled_error","axis_a_relative_error")
    for r in mesh:
        h=r["history"][-1]
        lines.append(f"| {r['configuration']['radial_cells']} | "+" | ".join(f"{h[k]:.7g}" for k in keys)+" |")
    lines += ["", "The Q1-gradient B error is approximately first order; flux, current and local heat refine approximately at second order. Axis-field error has its own slower convergence. Boundary current is the material-side trace at the sphere vertices, normalized by the analytic equatorial current magnitude, not divided by a vanishing axis current. Local heat compares exact cell-integrated heat on actual polygonal cells. Global energy/heat also retain geometric error relative to the true sphere.","",
              "| Case | B relative L2 | Current relative L2 | Local heat relative L2 | Accepted diffusion substeps / rejected trials | Response residual [J] |",
              "|---|---:|---:|---:|---:|---:|"]
    for r in decay["cases"]:
        h=r["history"][-1]
        lines.append(f"| {r['name']} | {h['field_relative_volume_L2']:.7g} | {h['current_relative_volume_L2']:.7g} | {h['local_heat_relative_volume_L2']:.7g} | {r['accepted_substeps']} / {r['rejected_trials']} | {r['response_residual_j']:.7g} |")
    lines += ["", "Time-half and tighter-tolerance controls use nr=32; matching, mode and quadrature controls use nr=16 with unchanged exterior radial spacing. `mode2` is the independently known second radial eigenmode, not numerically formed reversal. Its larger Q1-gradient field error remains visible even though current and heating errors are smaller. All configurations and entire sampled histories are archived.","",
              "## Stiff time integration: energy balance is insufficient","",
              "An independently diagonalized lowest spatial eigenmode gives an exact semidiscrete decay exp(-4). This isolates time integration from continuum spatial error. A single raw midpoint step returns amplitude -1/3: it falsely reverses the field while satisfying magnetic-loss/heat bookkeeping. The adaptive method retains positive decay and improves with tolerance:","",
              "| Local tolerance | Final amplitude / initial | Absolute amplitude error | Accepted substeps | Rejected trials |",
              "|---:|---:|---:|---:|---:|"]
    stiff=decay["stiff_control"]
    for r in stiff["adaptive"]:
        lines.append(f"| {r['tolerance']:g} | {r['amplitude']:.10g} | {r['amplitude_absolute_error']:.7g} | {r['accepted_substeps']} | {r['rejected_trials']} |")
    lines += ["", f"Exact amplitude is {stiff['exact_amplitude']:.10g}; the raw-step balance residual is only {stiff['raw_balance_residual_j']:.4g} J despite its incorrect sign. The controller accepts the uncorrected two-half-step solution and its physical heat. It does not add damping, clip the field or force an energy residual to zero. A local tolerance is not a global error or sign guarantee for arbitrarily small fields; topology claims need resolution and uncertainty controls.","",
              "## Actual moving resistive plasma","",
              "Each moving case begins in the analytic spherical FRC corresponding to its initial bias: a=-3B(0)(1-s²)/4 and p=(45/8)B(0)² r²(1-s²)/mu0. Positive eta then changes flux and deposits heat, so subsequent motion is resistive relaxation rather than pure equilibrium error. Source cases also apply B(t)=sqrt(mu0) exp(-(t-0.04)²/(2*0.04²)) T, ending at its peak at 0.04 s. Scalar D is prescribed mathematical physics, not a validated argon transport model.","",
              "| Case | D [m²/s] | Max speed [m/s] | Physical Ohmic heat [J] | Source work [J] | Response residual [J] | Outer steps | Diffusion substeps / rejects |",
              "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for name in [f"{kind}_{n}" for kind in ("frc","source") for n in (4,8,16,32)]+["source_wide_16","source_wide_32"]:
        r=rows[name];d=r["final"]
        lines.append(f"| {name} | {r['configuration']['diffusivity']:g} | {d['maximum_speed_m_s']:.7g} | {d['ohmic_heat_j']:.7g} | {d['source_work_j']:.8g} | {d['response_energy_residual_j']:.7g} | {d['steps']} | {d['diffusion_accepted_substeps']} / {d['diffusion_rejected_trials']} |")
    lines += ["", "Changing D from 0.1 to 1 m²/s changes physical diffusion/heating, and is not a numerical stabilization choice. The associated scale sqrt(D*0.04) is about 0.0632 or 0.2 m, respectively; this is only a diffusion-length estimate for a changing, moving geometry. The model has no exterior mass or Ohmic heat. Gas energy, field energy, physical heating and source work are separately recorded.","",
              "### Material-coordinate spatial comparisons","",
              "Corresponding material vertices supply velocity errors, using coarser inertial-mass weights. Actual Q1 fields/current are sampled at common reference-sector points with composite quadrature over every fine subcell. Density/pressure restrict by actual cell volume and integrated heat sums before division by volume. The two evolved coordinate maps can differ: these are material-coordinate state differences, not Cartesian pointwise errors. Changes subtract each mesh's initial state; total current is additionally compared without subtraction. All values below are maximum differences across the midpoint/final saved states, relative to the finer reference.","",
              "| Coarse / fine | Velocity [%] | Magnetic change [%] | Total current [%] | Current change [%] | Density change [%] | Pressure change [%] | Local heat density [%] |",
              "|---|---:|---:|---:|---:|---:|---:|---:|"]
    normkeys=("velocity_relative_mass_L2","magnetic_change_relative_volume_L2","current_relative_volume_L2",
              "current_change_relative_volume_L2","density_change_relative_volume_L2","pressure_change_relative_volume_L2",
              "local_heat_density_relative_volume_L2")
    mesh_names={f"{kind}_{n}" for kind in ("frc","source") for n in (4,8,16)}|{"source_wide_16"}
    for r in controls:
        if r["test"] in mesh_names and r["reference"].endswith(("_8","_16","_32")):
            lines.append(f"| {r['test']} / {r['reference']} | "+" | ".join(f"{100*max(t[k] for t in r['samples']):.6g}" for k in normkeys)+" |")
    source=next(r for r in controls if r["test"]=="source_16" and r["reference"]=="source_32")
    wide=next(r for r in controls if r["test"]=="source_wide_16")
    lines += ["", f"At D=0.1, the 16-to-32 velocity change is {100*max(t[normkeys[0]] for t in source['samples']):.5g}%; at D=1 it is {100*max(t[normkeys[0]] for t in wide['samples']):.5g}%. These adjacent-grid differences are not absolute errors or proof of a fully converged formation experiment. Current/thermal/density changes must be judged alongside velocity and the fixed-sphere analytic checks.","",
              "### Time, tolerance and magnetic-boundary controls","",
              "The following use the same physical model and matching material resolution. A lower outer timestep does not replace diffusion error control, and a tighter diffusion tolerance does not replace material-wave/source-time resolution. Matching radius is varied at fixed exterior radial spacing.","",
              "| Test / reference | Velocity relative L2 [%] | Velocity RMS difference [m/s] | Magnetic-change L2 [%] | Total-current L2 [%] | Local-heat-density L2 [%] |",
              "|---|---:|---:|---:|---:|---:|"]
    for r in controls:
        if r["test"] in mesh_names and r["reference"].endswith(("_8","_16","_32")):
            continue
        samples=r["samples"]
        values=[100*max(t[normkeys[0]] for t in samples),max(t["velocity_rms_difference_m_s"] for t in samples),
                100*max(t[normkeys[1]] for t in samples),100*max(t[normkeys[2]] for t in samples),100*max(t[normkeys[-1]] for t in samples)]
        lines.append(f"| {r['test']} / {r['reference']} | "+" | ".join(f"{v:.7g}" for v in values)+" |")
    lines += ["", "| Moving fixture | Outer time cap [s] | Steps | Response residual [J] | Ohmic heat [J] |",
              "|---|---:|---:|---:|---:|"]
    for kind in ("frc","source"):
        for name in (f"{kind}_time_0.001",f"{kind}_8",f"{kind}_time_0.00025"):
            r=rows[name];d=r["final"]
            lines.append(f"| {kind}, nr=8 | {r['configuration']['max_timestep']:g} | {d['steps']} | {d['response_energy_residual_j']:.8g} | {d['ohmic_heat_j']:.8g} |")
    mass=max(abs(h["mass_kg"]/r["initial"]["mass_kg"]-1) for r in rows.values() for h in r["history"])
    all_history=[h for r in rows.values() for h in r["history"]]
    lines += ["", "## Interpretation and remaining requirements","",
              f"All sampled material mass changes are {mass:.4g} relative. Maximum sampled resistive-substep energy residual is {max(h['maximum_resistive_substep_balance_residual_j'] for h in all_history):.4g} J, while the full moving response has its separate finite-time error above. Maximum sampled vacuum relative residual is {max(h['vacuum_relative_residual'] for h in all_history):.4g}. Minimum sampled density and pressure are {min(h['minimum_density_kg_m3'] for h in all_history):.7g} kg/m³ and {min(h['minimum_pressure_pa'] for h in all_history):.7g} Pa. These sampled extrema are not bounds between outputs. Heat is nonnegative and actual cumulative cell heat is archived; there is no pressure floor.","",
              "Independent polynomial/volume checks validate the material-field sampler and conservative cell restrictions used in these comparisons. Atomic-failure tests cover a rejected second complete resistive half-step, after trial heating and material motion, and diffusion trial exhaustion. The accepted state, labels, entropy, heat, time and source work remain unchanged on failure.","",
              "This extends real moving-interface dynamics to physical resistive induction and local heating. It does not establish physical transport validity, surface thermodynamic closure, arbitrary coil geometry, shocks/remapping, long-time mesh robustness, Hall effects, ionization, particle confinement or FRC formation. The initial reversed-field topology is prescribed. Numerical field accuracy, current-layer resolution, actual boundary/thermal behavior and model validity remain independent requirements.","",
              "Next add spatially varying prescribed circular coils with the full geometric source-force derivative and source-work verification, while addressing the measured field/current/flow error and physical transport/Hall applicability. The dynamic experiment interface must preserve those validity limits rather than unlocking an unverified formation claim.","",
              "Full configurations, accepted histories, evolving labels/entropy, field/current/heat arrays, environment and hashes are in `outputs/resistive_moving_verification`. Reproduce with `python -m examples.resistive_decay_measurements`, `python -m examples.resistive_dynamic_measurements --reuse`, then `python -m examples.verify_resistive_moving`. The finest moving cases take minutes. Report/analysis/render sources are separately hashed from the frozen scientific measurements.","",
              "![Resistive accuracy](../outputs/resistive_moving_verification/resistive_verification.png)","",
              "![Computed resistive fields](../outputs/resistive_moving_verification/resistive_fields.png)"]
    Path("docs/RESISTIVE_MOVING_VERIFICATION.md").write_text("\n".join(lines)+"\n")


def verify(directory="outputs/resistive_moving_verification"):
    directory=Path(directory);decay,rows,arrays=read_study(directory);controls=comparisons(rows,arrays)
    files=[Path(f"examples/{name}.py") for name in ("resistive_analysis","resistive_plot","verify_resistive_moving",
            "verification_records","moving_plot","moving_analysis")]
    result={"completed":True,"provenance":decay["provenance"],"report_environment":{"matplotlib":matplotlib.__version__},
            "report_sources":{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            "decay":decay,"dynamic":list(rows.values()),"comparisons":controls}
    write_record(directory/"metrics.json",result)
    report(decay,rows,controls);error_figure(directory,decay,rows);field_figure(directory,rows,arrays)
    print("Verified complete resistive source/configuration/field hashes; report and figures written",flush=True)


if __name__=="__main__":
    verify()
