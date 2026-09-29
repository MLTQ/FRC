"""Assemble verified moving-material records, report and scientific figures."""
from pathlib import Path
import hashlib
import matplotlib
from .moving_analysis import read_study,analyze,nonuniform_advection
from .moving_plot import plot_errors,plot_fields
from .verification_records import write_record


def report(static,cases,analysis,advection):
    lines=["# Moving material plasma / infinite vacuum verification", "",
           "This study verifies an experimental ideal, adiabatic material mesh with an actual moving plasma surface and no exterior fluid. A weighted regular-axis magnetic finite-element solve is matched continuously to analytic decaying modes at infinity. The [method derivation](MOVING_PLASMA.md) states equations, quadrature, source work, advection approximations and limitations. Existing fixed-chamber methods and defaults are unchanged; Stage 2 formation remains gated.", "",
           "The 26 evolution cases and 27 static evaluations use SI unit-scale mathematical fixtures, not a selected gas/temperature/transport regime. Every evolution begins independently. All complete configurations, accepted histories, source/environment hashes and lossless field archives are retained in `outputs/moving_verification`. Old fixed-boundary archives retain their own historical sources.", "",
           "## Independent magnetic energy and generalized force", "",
           "For a unit initial sphere with Bp0=sqrt(mu0), uniform bias=0.3 Bp0 and homologous material scale lambda, fixed flux gives Bp=Bp0/lambda². Integrating the interior uniform induced field and the exterior dipole gives E=pi(1/lambda-0.6 lambda+0.09 lambda³) J. Its independent derivative is pi(-1/lambda²-0.6+0.27 lambda²). No assembled operator is used to construct these targets. This electromagnetic fixture is not a valid free-gas equilibrium with its nonzero normal field and tangential-field jump.", "",
           "| Scale | Radial cells | Energy relative error | Generalized-force relative error |",
           "|---|---:|---:|---:|"]
    for r in static["homologous"]:
        lines.append(f"| {r['scale']:.1f} | {r['radial_cells']} | {r['energy_relative_error']:.7g} | {r['denergy_dscale_relative_error']:.7g} |")
    lines += ["", "The finite material sphere is polygonal and has second-order geometric error; the distant matching sphere is exactly curved. Force is the analytic derivative of the full discrete potential, including vacuum mesh geometry, material flux coefficients and the regular-axis chain rule. An arbitrary two-dimensional deformation is separately checked against a symmetric perturbation of that complete potential. Agreement of two discrete energy quadratures does not remove continuum error.", "",
              "### Translation, modes and matching", "",
              "Axial translation in uniform bias leaves exact full-space energy unchanged. The finite magnetic mesh breaks this symmetry slightly; its net axial force is a numerical error, not a physical centering force.", "",
              "| Radial cells | Shift [m] | Spurious axial force [N] | Energy relative error |",
              "|---:|---:|---:|---:|"]
    for r in static["translation"]:
        lines.append(f"| {r['radial_cells']} | {r['axial_shift_m']} | {r['spurious_axial_force_n']:.7g} | {r['energy_relative_error']:.7g} |")
    lines += ["", "The next controls use 12 radial cells and a 0.4 m axial translation. Moving the matching sphere preserves exterior radial spacing, but changes exterior element shapes; matching sensitivity can persist even with a tiny elliptic residual.", "",
              "| Control | Value | Energy relative error | Generalized-force relative error | Axial force [N] |",
              "|---|---:|---:|---:|---:|"]
    for key,selector in (("modes","modes"),("matching","matching_radius_m"),("quadrature","quadrature")):
        for r in static[key]:
            lines.append(f"| {key} | {r[selector]:g} | {r['energy_relative_error']:.7g} | {r['denergy_dscale_relative_error']:.7g} | {r['spurious_axial_force_n']:.7g} |")
    maximum_vacuum=max(r["vacuum_relative_residual"] for key in ("homologous","translation","modes","matching","quadrature") for r in static[key])
    lines += ["", f"Maximum stationary-vacuum relative residual in the static study is {maximum_vacuum:.4g}. Increasing modes eventually reaches spatial-error limits. Exact magnetic axial-momentum conservation is not claimed.", "",
              "## Actual free gas expansion", "",
              "For gamma=5/3, rho0=1 kg/m³, p0=1-s² Pa, B=0 and zero initial velocity inside the unit sphere, lambda''=2/lambda³ gives lambda=sqrt(1+2t²), x=lambda X, v=2t X/lambda. Surface pressure is zero and the exterior has no mass. The table compares actual evolution to this independent solution at 0.1 s. Velocity uses positive initial nodal-mass weights; pressure compares exact cell averages on the actual final geometry.", "",
              "| Radial cells | Steps | Velocity relative L2 | Maximum surface-position error [m] | Density relative L2 | Pressure cell-average relative L2 | Energy residual [J] |",
              "|---:|---:|---:|---:|---:|---:|---:|"]
    for n in (4,8,16,32):
        r=cases[f"expansion_{n}"];e=r["final_errors"];d=r["final"]
        lines.append(f"| {n} | {d['steps']} | {e['velocity_relative_mass_L2']:.7g} | {e['maximum_boundary_position_error_m']:.7g} | {e['density_relative_volume_L2']:.7g} | {analysis['expansion_pressure_cell_average_error'][r['name']]:.7g} | {d['response_energy_residual_j']:.7g} |")
    lines += ["", "The exact final radius is 1.009950494 m. Initial polygonal volume error and P0 entropy projection are retained; no correction makes the reference exact on a coarse mesh. Refinement of the velocity and surface errors provides stronger evidence than exact material mass alone.", "",
              "## Initialized spherical FRC evolution", "",
              "Initialize the independently derived [spherical FRC](FREE_BOUNDARY_REFERENCE.md) with B0=sqrt(mu0), unit radius, rho=1 kg/m³ and gamma=5/3. Inside, a=-3 B0(1-s²)/4 and p=(45/8)r²(1-s²) Pa. Outside, uniform bias plus the solved decaying field supplies vacuum traction. The continuum state has zero acceleration; all computed motion below is numerical drift. Its exact internal and induced energies are 9pi/7 and 10pi/7 J.", "",
              "| Radial cells | Initial mass-RMS acceleration [m/s²] | Max speed at 0.04 s [m/s] | Max surface displacement [m] | Energy residual [J] |",
              "|---:|---:|---:|---:|---:|"]
    for n in (4,8,16,32):
        r=cases[f"frc_{n}"];e=r["final_errors"];d=r["final"]
        lines.append(f"| {n} | {r['initial_errors']['initial_or_current_acceleration_mass_rms_m_s2']:.7g} | {d['maximum_speed_m_s']:.7g} | {e['maximum_boundary_displacement_m']:.7g} | {d['response_energy_residual_j']:.7g} |")
    lines += ["", "This is evolution of an initialized equilibrium, not creation of an FRC. The finite-time drift and mesh/time controls do not establish long-time stability. Density is a mathematical choice, not a validated thermodynamic closure near the zero-pressure axis and surface.", "",
              "## Time and prescribed-bias work", "",
              "At nr=8, time-cap refinement uses identical spatial data. Material CFL can impose a smaller step; accepted step counts below show the actual amount of work. Every 0.004 s diagnostic boundary is reached exactly. The response residual is kinetic plus gas internal plus all-space induced magnetic energy, minus initial energy and integrated prescribed-source work; it is diagnosed without repair.", "",
              "| Fixture | Time cap [s] | Steps | Source work [J] | Response residual [J] | Maximum speed [m/s] |",
              "|---|---:|---:|---:|---:|---:|"]
    for kind in ("expansion","frc","source"):
        for dt in (.002,.001,.0005,.00025):
            r=cases[f"{kind}_time_{dt}"];d=r["final"]
            lines.append(f"| {kind} | {dt:g} | {d['steps']} | {d['source_work_j']:.8g} | {d['response_energy_residual_j']:.8g} | {d['maximum_speed_m_s']:.7g} |")
    lines += ["", "The source fixture starts in the FRC equilibrium appropriate to B(0), then applies B(t)=sqrt(mu0) exp(-(t-0.04)²/(2*0.04²)) T through its peak at 0.04 s. The source is a prescribed spatially uniform magnetic bias, not arbitrary circular coils or a circuit model. The exact source-power formula is checked independently by differentiating magnetic energy at fixed material coordinates. Its driven flow has no analytic reference, so adjacent-grid differences remain necessary.", "",
              "| Source radial cells | Max speed [m/s] | Max material displacement [m] | Source work [J] | Response residual [J] |",
              "|---:|---:|---:|---:|---:|"]
    for n in (4,8,16):
        r=cases[f"source_{n}"];d=r["final"]
        lines.append(f"| {n} | {d['maximum_speed_m_s']:.7g} | {r['final_errors']['maximum_material_displacement_m']:.7g} | {d['source_work_j']:.8g} | {d['response_energy_residual_j']:.7g} |")
    lines += ["", "### Saved-state comparisons", "",
              "Differences use corresponding material nodes at matching physical saved times. The table reports the maximum across midpoint and final states, with the first row's mass weights; relative differences in tiny artificial equilibrium motion are not relative physical error.", "",
              "| Test / reference | Velocity relative L2 | Velocity RMS difference [m/s] | Maximum position difference [m] |",
              "|---|---:|---:|---:|"]
    for r in analysis["comparisons"]:
        values=[max(t[key] for t in r["samples"]) for key in ("velocity_relative_mass_L2","velocity_mass_rms_difference_m_s","maximum_position_difference_m")]
        lines.append(f"| {r['test']} / {r['reference']} | "+" | ".join(f"{v:.7g}" for v in values)+" |")
    driven=next(r for r in analysis["comparisons"] if r["test"]=="source_8" and r["reference"]=="source_16")
    flow_change=100*max(t["velocity_relative_mass_L2"] for t in driven["samples"])
    lines += ["",f"The driven 8-to-16 mesh change is {flow_change:.5g}% in velocity L2. This is not a converged driven-flow result, even when source-work error is much smaller. Ideal surface-current response launches a narrow moving disturbance; spatial refinement and the later independently verified resistive extension must resolve it before interpreting formation. The current coarse-grid drive is a source/work and sensitivity fixture only."]
    lines += ["", "## Flux inside material cells", "",
              "The solver preserves psi at each moving non-axis vertex. It does not preserve every Q1-interpolated material-point psi exactly. For the independent prescribed deformation r_new=r(1+0.2z), z_new=z+0.1r² and initial a=0.5(1-0.2r²+0.1z) T, the actual interior discrepancy is:", "",
              "| Radial cells | Material-point flux relative L2 drift |","|---:|---:|"]
    for r in advection["rows"]:
        lines.append(f"| {r['radial_cells']} | {r['material_point_psi_relative_L2']:.8g} |")
    lines += ["", "This second-order spatial defect is numerical interpolation error, not physical resistivity or reconnection. For the evolved FRC/source cases:", "",
              "| Case | Max vertex psi drift [Wb/radian] | Max material-point flux relative L2 drift |","|---|---:|---:|"]
    for name in [f"frc_{n}" for n in (4,8,16,32)]+[f"source_{n}" for n in (4,8,16)]:
        vertex=cases[name]["final_errors"]["nodal_psi_max_drift_wb_per_radian"]
        interior=max(r["material_point_psi_relative_L2"] for r in analysis["interpolated_material_flux"][name])
        lines.append(f"| {name} | {vertex:.5g} | {interior:.7g} |")
    histories=[h for r in cases.values() for h in r["history"]]
    mass_error=max(abs(h["mass_kg"]/r["initial"]["mass_kg"]-1) for r in cases.values() for h in r["history"])
    lines += ["", "## Budgets, interpretation and remaining work", "",
              f"All sampled material mass changes are {mass_error:.4g} relative. The maximum sampled vacuum relative residual is {max(h['vacuum_relative_residual'] for h in histories):.4g}. Minimum sampled density and pressure are {min(h['minimum_density_kg_m3'] for h in histories):.7g} kg/m³ and {min(h['minimum_pressure_pa'] for h in histories):.7g} Pa, without floors. These sampled extrema do not assert bounds between outputs. Every accepted step validates positive geometry and finite state; rejected proposals retain the prior accepted state, time and work.", "",
              "The evolving material interface is a new capability. The measured ideal interpolation defect, discrete equilibrium drift, spurious translation force and matching dependence remain explicit. Small budget and elliptic residuals do not prove complete physical accuracy. Q1/P0 material meshes can distort and have low-energy modes; there is no shock capture or remapping in this reference.", "",
              "Next extend the compatible ideal force/vacuum treatment to physical resistive induction and positive Ohmic heating with independent diffusion/energy tests, then spatial coil sources and resolved driven-flow comparisons. Continue transport/ionization/Hall applicability, actual geometry/end conditions, shock/mesh robustness and converged formation diagnostics. An initialized reversed-field torus, a closed contour, numerical survival and particle confinement are distinct claims.", "",
              "Reproduce measurements with `python -m examples.moving_static_measurements` and `python -m examples.moving_measurements --reuse`; assemble only complete source-matched records with `python -m examples.verify_moving`. The finest evolution cases can take several minutes. Sources for continuum concepts and the complete derivation are linked in MOVING_PLASMA.md; no published stability guarantee is applied to this discretization.", "",
              "![Independent convergence](../outputs/moving_verification/convergence.png)", "",
              "![Computed material fields](../outputs/moving_verification/material_fields.png)"]
    Path("docs/MOVING_PLASMA_VERIFICATION.md").write_text("\n".join(lines)+"\n")


def verify(directory="outputs/moving_verification"):
    directory=Path(directory);static,cases,arrays=read_study(directory)
    analysis=analyze(cases,arrays);advection=nonuniform_advection()
    files=[Path(f"examples/{name}.py") for name in ("moving_analysis","moving_plot","verify_moving","verification_records")]
    m={"provenance":static["provenance"],"report_environment":{"matplotlib":matplotlib.__version__},
       "report_sources":{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
       "completed":True,"static":static,"cases":list(cases.values()),"analysis":analysis,"nonuniform_advection":advection}
    write_record(directory/"metrics.json",m)
    report(static,cases,analysis,advection)
    plot_errors(directory,static,cases);plot_fields(directory,cases,arrays)
    print("Verified all moving-plasma source/configuration/field hashes; report and figures written",flush=True)


if __name__=="__main__":
    verify()
