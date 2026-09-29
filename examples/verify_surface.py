"""Assemble complete surface studies without relabeling mesh differences as errors."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
from .surface_analysis import read_study,comparison_pairs
from .surface_comparisons import measurements
from .surface_forces import forces,initial_compatibility
from .surface_plot import profiles,accuracy
from .verification_records import write_record


def historical_baseline(arrays):
    result=[];root=Path("outputs/resistive_moving_verification/dynamic")
    for n in (8,16,32):
        old=json.loads((root/f"source_{n}.json").read_text());path=root/old["fields"]["file"]
        if old["provenance"]["source_sha256"]!="452e51b4405b8a44fc6191f6716fa4e19e5798ba367bb84a0953ce3c6703503b":
            raise ValueError("Unexpected historical baseline source")
        if hashlib.sha256(path.read_bytes()).hexdigest()!=old["fields"]["sha256"]:
            raise ValueError("Changed historical field archive")
        with np.load(path,allow_pickle=False) as data:
            current=arrays[f"source_g1_{n}"]
            changed=[key for key in data.files if not np.array_equal(data[key],current[key])]
            if set(data.files)!=set(current):raise ValueError("Baseline field contract changed")
        result.append({"radial_cells":n,"historical_source_sha256":old["provenance"]["source_sha256"],
                       "changed_arrays":changed,"bitwise_identical":not changed})
    return result


def report(directory,result):
    rows={r["name"]:r for r in result["dynamic"]};refs=result["references"]["cases"]
    lines=["# Resolving the moving resistive surface layer","",
      "This study adds explicit radial grading to the actual material/vacuum mesh while retaining the same physical D=eta/mu0=0.1 m²/s. Material shells use R[1-(1-xi)^p], and vacuum offsets use (Rmatch-R)xi^p. The default p=1 retains uniform spacing. This is an initial mesh choice, not adaptive remapping, artificial resistance or a transport closure.","",
      "All 21 moving cases and 13 independent component cases are complete with source/configuration/field hashes checked. The extra finest-mesh time control was added after the first control exposed time sensitivity, using the same frozen measurement entry point and fully recorded typed inputs; its reproduction helper is separately hashed. The cases use the preceding mathematical FRC and uniform-bias fixtures, rho=1 kg/m³ and gamma=5/3. Reversed flux is initialized. These are not demonstrated formation experiments or validated species parameters.","",
      "## Independent graded-mesh checks","",
      "Physical volume and second moments match independent disk integration of the revolved polygon. Actual free-gas expansion is checked against lambda=sqrt(1+2t²), with p0=1-s² and rho0=1. Fixed conducting-sphere decay independently checks full profiles, current, local Ohmic heat and surface current. The gas test has no heat source; the diffusion test has fixed geometry. Neither alone validates a heated free surface.","",
      "| Decay case | B relative L2 | Current relative L2 | Local heat relative L2 | Surface-current scaled maximum error |",
      "|---|---:|---:|---:|---:|"]
    for r in refs:
        if r["name"].startswith("decay"):
            h=r["history"][-1];keys=("field_relative_volume_L2","current_relative_volume_L2","local_heat_relative_volume_L2","maximum_boundary_current_scaled_error")
            lines.append(f"| {r['name']} | "+" | ".join(f"{h[k]:.7g}" for k in keys)+" |")
    lines += ["", "| Free-gas case, t=0.04 s | Velocity relative mass L2 | Maximum surface-position error [m] | Response-energy residual [J] |",
              "|---|---:|---:|---:|"]
    for r in refs:
        if r["name"].startswith("expansion"):
            e=r["errors"];lines.append(f"| {r['name']} | {e['velocity_relative_mass_L2']:.7g} | {e['maximum_boundary_position_error_m']:.7g} | {r['final']['response_energy_residual_j']:.7g} |")
    lines += ["", "Current and local heat retain approximately second-order refinement, while Q1-gradient B remains about first order. Grading trades interior resolution for surface resolution; the smooth gas-expansion error therefore also needs measurement. No automatic accuracy claim follows from a smaller last cell.","",
      "## Actual moving fields and surface motion","",
      "All moving runs end at 0.04 s. Source cases use B(t)=sqrt(mu0) exp(-(t-0.04)²/(2*0.04²)) T; unforced FRC cases hold B=sqrt(mu0) T. Positive resistivity changes flux and deposits heat. The initial pressure vanishes at the analytic surface while current does not, so surface thermal evolution is a distinct physical and numerical requirement.","",
      "| Case | nr × nt | Grading | Peak speed [m/s] | Equatorial outward surface speed [m/s] | Minimum density [kg/m³] | Physical heat [J] | Steps |",
      "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows.values():
        c=r["configuration"];d=r["final"]
        lines.append(f"| {r['name']} | {c['radial_cells']} × {c['angular_cells']} | {c['radial_grading']:g} | {d['maximum_speed_m_s']:.7g} | {result['surface_speed'][r['name']]:.7g} | {d['minimum_density_kg_m3']:.7g} | {d['ohmic_heat_j']:.7g} | {d['steps']} |")
    lines += ["", "The surface and the immediately interior material can move in opposite directions. The current, density, pressure and complete surface profiles are archived and plotted; bulk flow norms must not conceal a poorly resolved boundary. No floor, exterior gas or energy repair is introduced.","",
      "## Comparison on common initial physical material coordinates","",
      "The sampler inverts each initial straight-sector map, including its graded radial shells, and evaluates the actual Q1 fields at the same initial (r,z). Velocity uses initial physical mass weights; magnetic/current/thermal comparisons use evolved test-volume weights. Field, density and pressure changes subtract each mesh's own initial sample. This differs from the prior coarse-vertex velocity norm, so the percentages below must not replace or be equated to that report's 63.821% number. Thermal fields also differ in definition: the earlier report conservatively restricted fine-cell values to coarse cells; this pointwise P0 comparison retains differences inside those coarse cells.","",
      "The angular-coarse polygon lies inside the angular-fine polygon. Its support omits a reported fraction of the finer initial volume. A separate maximum boundary-vertex velocity difference compares shared initial surface locations and retains the surface-motion uncertainty. Differences from a finer run are not absolute continuum errors.","",
      "| Test / reference | Velocity mass L2 [%] | Magnetic change [%] | Total current [%] | Pressure change [%] | Density change [%] | Local heat density [%] | Max boundary-speed difference [m/s] | Omitted initial volume [%] |",
      "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    keys=("velocity_relative_mass_L2","magnetic_change_relative_volume_L2","total_current_relative_volume_L2","pressure_change_relative_volume_L2","density_change_relative_volume_L2","local_heat_density_relative_volume_L2")
    for c in result["comparisons"]:
        values=[100*max(s[k] for s in c["samples"]) for k in keys]
        boundary=max(s["maximum_boundary_velocity_difference_m_s"] for s in c["samples"])
        lines.append(f"| {c['test']} / {c['reference']} | "+" | ".join(f"{v:.6g}" for v in values)+f" | {boundary:.7g} | {100*c['omitted_reference_initial_volume_fraction']:.6g} |")
    lines += ["", "Each table entry is the maximum across the midpoint/final saved times. Pressure/density/heat remain actual discontinuous cell fields; no smoothing is used. Initial material-cell intersections form a common integration mesh. Positive triangle rules integrate inside each overlap, and actual deformed Jacobian/radius ratios provide volume weights. The integrated overlap volume is checked against the whole test polygon. Comparison quadrature is checked independently by increasing its one-dimensional Gauss order from 4 to 6:","",
              "| Norm | Maximum relative change in the measured difference [%] |","|---|---:|"]
    for key,value in result["comparison_quadrature_relative_changes"].items():
        lines.append(f"| {key} | {100*value:.6g} |")
    lines += ["", "These integration sensitivities are part of the diagnostic uncertainty, especially for piecewise constant thermal fields. They are not time-step errors. Explicit time/CFL, angular, diffusion tolerance, magnetic-matching and physical quadrature controls are distinct rows above.","",
      "The finest grade-2 mesh still changes by 15.8362% in velocity when angular resolution doubles. Halving CFL and the outer step cap gives 5.39874%, then 2.75800% velocity differences on the same mesh. Their ratio is about 1.96, not evidence of second-order asymptotic convergence; current differences do not decrease monotonically. These controls retain their fixed diffusion tolerance. Surface-speed differences remain 0.04457785 then 0.01810780 m/s. No time-extrapolated solution or resolved boundary is claimed.","",
      "## Force cancellation and accounting","",
      "The following final equatorial surface accelerations come from restoring the saved accepted state. The recovered-volume J cross B column separately projects weak current times Q1 B. Both magnetic columns are discrete constructions; their agreement is a diagnostic, not an independent continuum-force proof.","",
      "| Case | Pressure [m/s²] | Magnetic virtual work [m/s²] | Recovered volume J cross B [m/s²] | Net [m/s²] |",
      "|---|---:|---:|---:|---:|"]
    for name,f in result["forces"].items():
        keysf=("pressure","magnetic","weak_current","net")
        lines.append(f"| {name} | "+" | ".join(f"{f[k+'_radial_acceleration_m_s2'][-1]:.7g}" for k in keysf)+" |")
    history=[h for r in rows.values() for h in r["history"]]
    mass=max(abs(h["mass_kg"]/r["initial"]["mass_kg"]-1) for r in rows.values() for h in r["history"])
    lines += ["",f"Maximum sampled relative mass change is {mass:.4g}; resistive-stage energy residual is {max(h['maximum_resistive_substep_balance_residual_j'] for h in history):.4g} J; vacuum relative residual is {max(h['vacuum_relative_residual'] for h in history):.4g}. Whole moving-response time error is independently retained in every history and is not repaired. Physical heat is nonnegative in all saved cells. These sampled statements are not bounds between outputs.","",
      f"The new uniform nr=8/16/32 runs have {'bitwise identical' if all(r['bitwise_identical'] for r in result['historical_baseline']) else 'changed (see metrics)'} field/state arrays compared with the prior verified uniform-grid checkpoint. Historical results keep their historical fingerprints.","",
      "## Remaining scope","",
      "A separate initial-time compatibility audit derives a positive equatorial Joule pressure source of 3.75 Pa/s (unforced FRC) and 1.379547904 Pa/s (initial source fixture). If B is continuous, no additional surface stress/heat term is present and p=0 remains a smooth material trace with bounded velocity gradient, its material time derivative would instead be zero. Thus these initial heated-edge data fail that conditional smooth-boundary compatibility requirement. Independent one-sided Maxwell derivatives check the current used in the calculation. See [the derivation, assumptions and primary sources](SURFACE_COMPATIBILITY.md). This does not prove nonexistence of every weak solution or excuse a numerical error; it changes how these initialized resistive runs can be used as convergence targets.","",
      "Surface resolution is now an explicit, independently checked experiment control. The measured changes and boundary-speed differences remain the evidence of what is and is not resolved. This does not establish heated-surface continuum closure, shock/remapping robustness, temperature-dependent transport, ionization or Hall validity, actual chamber/end boundaries, particle confinement or FRC formation. Spatial coil coupling remains the next source capability, with its independent force reference in MOVING_COIL_PLAN.md. Stage 2 remains gated.","",
      "Reproduce the component records with `python -m examples.surface_measurements --references`, the original 20 moving records with `python -m examples.surface_measurements --reuse`, the additional time control with `python -m examples.surface_time_control`, and this assembly with `python -m examples.verify_surface`. The finest angular/time controls take many minutes. All source/configuration/field hashes and separate analysis/render/reproduction sources are retained in `outputs/surface_verification/metrics.json`.","",
      "![Computed surface profiles](../outputs/surface_verification/surface_profiles.png)","",
      "![Surface accuracy](../outputs/surface_verification/surface_accuracy.png)"]
    Path("docs/SURFACE_VERIFICATION.md").write_text("\n".join(lines)+"\n")


def verify(directory="outputs/surface_verification"):
    directory=Path(directory);refs,rows,arrays=read_study(directory)
    measured=measurements(directory,rows,arrays)
    if len(measured)!=len(comparison_pairs()):raise ValueError("Missing declared comparison")
    controls=[r["order4"] for r in measured]
    qchanges={}
    for item in measured:
        entry,fine=item["order4"],item["order6"]
        for low,high in zip(entry["samples"],fine["samples"]):
            for key in low:
                if "relative_" in key:
                    qchanges[key]=max(qchanges.get(key,0.),abs(low[key]-high[key])/max(abs(high[key]),1e-14))
    names=("source_g1_16","source_g1.5_32","source_g2_16","source_g2_32","source_g2_32_angular")
    force_rows={name:forces(name,arrays) for name in names};speeds={}
    for name,row in rows.items():
        a=arrays[name];node=a["t10_boundary"][row["configuration"]["angular_cells"]//2]
        speeds[name]=float(a["t10_velocity"][node,0])
    files=[Path(f"examples/{n}.py") for n in ("surface_analysis","surface_comparisons","surface_sampling","material_overlap","surface_forces","surface_plot","verify_surface","surface_time_control","spherical_frc_reference","verification_records")]
    result={"completed":True,"provenance":refs["provenance"],"report_environment":{"matplotlib":matplotlib.__version__},
            "report_sources":{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            "references":refs,"dynamic":list(rows.values()),"comparisons":controls,
            "comparison_quadrature_relative_changes":qchanges,"forces":force_rows,"surface_speed":speeds,
            "historical_baseline":historical_baseline(arrays),"initial_surface_compatibility":initial_compatibility()}
    profiles(directory,rows,arrays,force_rows);accuracy(directory,refs,controls);report(directory,result)
    write_record(directory/"metrics.json",result)
    print("Completed surface study checked; report and figures written",flush=True)


if __name__=="__main__":
    verify()
