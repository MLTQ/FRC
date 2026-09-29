"""Assemble source, skin, end and fine-shock evidence for boundary order reduction."""
from pathlib import Path
import hashlib
import json
from .boundary_analysis import read_case, compare, compare_ends, summarize
from .boundary_pulse_measurements import protocol, source_record
from .boundary_plot import plot
from .cylindrical_skin_reference import CylindricalSkin
from .verification_records import write_record


def report(m):
    lines = ["# Boundary accuracy verification", "",
             "This study tests a fixed first-order wall/end strip against the compatible second-order method. The shared conservative/CT assembly, physical resistivity, natural vacuum work and SSP RK3 are unchanged. The strip remains an explicit experimental probe; production defaults are unchanged. Fixed order reduction is recorded separately from adaptive fallback.", "",
             "## Original shock, unchanged initial state", "",
             "The independent central shock reference and the 0.02 s end time are identical to the compatible-energy study. Full-space response error remains distinct from exactly conservative fluid energy-flux accounting.", "",
             "| Method | Radial cells | B L1/B1 | Pressure L1 [Pa] | Response residual [J] | Min density [kg/m³] | Additional adaptive limited stages |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for method, rows in m["shocks"].items():
        for r in rows:
            d, e = r["diagnostics"], r["errors"]
            lines.append(f"| {method} | {r['radial_cells']} | {e['Bz_L1_over_upstream_B']:.7g} | {e['pressure_L1_pa']:.7g} | {d['response_energy_balance_residual_j']:.7g} | {d['minimum_density_kg_m3']:.6g} | {d['limited_euler_stages']} |")
    fine = m["shocks"]["strip"][-1]
    p = fine["boundary_order_probe"]
    lines += ["", f"The nr=384 strip has {p['forced_cell_count']} forced cells and {p['affected_cell_count']} cells sharing selected faces/corners. Its adaptive count excludes these fixed selections. The earlier nr=192 time control decreases response residual from 3.973689e-7 to 4.788792e-8 J on CFL halving. The new fine-grid result alone does not replace a separate nr=384 temporal study.", "",
              "## Independent stationary-cylinder skin test", "",
              "For an infinite cylinder, Bz solves B_t=D(B_rr+B_r/r), with prescribed exterior Bz=Re(Bs exp(iωt)), regular axis and continuous boundary B. The independent solution uses k=(1+i)/δ, δ=√(2D/ω), Bhat=Bs I0(kr)/I0(kR), psihat=Bs r I1(kr)/(k I0(kR)), and Jhat=-Bs k I1(kr)/(μ0 I0(kR)). Outside, J=0 and psi continues as psi(R)+Bs(r²-R²)/2. The exact periodic solution initializes every case; eight phases over one period are measured.", "",
              "The magnetic-diffusion and skin-effect basis follows [MIT, magnetic diffusion laws](https://web.mit.edu/6.013_book/www/chapter10/10.5.html) and [MIT, skin effect](https://www.mit.edu/course/6/6.013_book/www/chapter10/10.7.html). The cylindrical solution is derived here by separation. Independent finite differences check Maxwell/Ohm/diffusion identities, and Gauss volume integration of mean ηJ² agrees with analytic inward surface Poynting power.", "",
              "This is a kinematic conductor check: gas is held at rho=1 kg/m³, p=1 Pa and v=0 to supply the declared production LLF speed. Standard/strip cases use the actual production reconstruction, CT electric fields and wall numerical resistance. Only the invariant-z central stencil is evolved; end/corner dynamics and gas motion are not tested here. The physical-only control omits both UCT and numerical interface diffusion. It is not a full MHD option.", "",
              "| Method | Radial cells | δ / grid spacing | Maximum B L2 / Bs | Maximum wall J error / (Bs/μ0δ) | Maximum flux error / (Bs R²/2) |",
              "|---|---:|---:|---:|---:|---:|"]
    for r in m["skin"]:
        if r["cfl"] == .2:
            lines.append(f"| {r['method']} | {r['radial_cells']} | {r['cells_per_skin_depth']:.5g} | {r['maximum_B_error']:.7g} | {r['maximum_wall_J_error']:.7g} | {r['maximum_flux_error']:.7g} |")
    ref = CylindricalSkin()
    lines += ["", f"R=1 m, D=0.1 m²/s, f=10 Hz give δ={ref.skin_depth:.7g} m. B errors compare exact cylindrical cell averages; norms use fixed applied-field scales, not a near-zero instantaneous norm. Mean analytic physical dissipation is {ref.mean_ohmic_power_per_length():.9g} W/m for the recorded tiny field amplitude. Surface current has a first-order pointwise stencil and is separately measured. Partial cancellation can make a numerically diffusive result closer than physical-only on a particular grid; that is not additional physical transport.", "",
              "| Method, nr=64 | CFL 0.2 maximum B error | CFL 0.1 maximum B error | Relative change [%] |",
              "|---|---:|---:|---:|"]
    for method in ("physical", "standard", "strip"):
        a, b = [next(r for r in m["skin"] if r["method"] == method and r["radial_cells"] == 64 and r["cfl"] == c) for c in (.2, .1)]
        lines.append(f"| {method} | {a['maximum_B_error']:.8g} | {b['maximum_B_error']:.8g} | {100*(a['maximum_B_error']/b['maximum_B_error']-1):.5g} |")
    lines += ["", "## Identical finite-cylinder coil pulses", "",
              "Both methods use the original two-coil pulse, fixed gas chamber, infinite multipole magnetic exterior, physical D=0.1 m²/s and 0.04 s rise. Source amplitude, density, pressure, and all constitutive/numerical controls are recorded. Full conserved, primitive, CT and current fields are saved at 0.04, 0.06 and 0.10 s in lossless hashed NPZ files.", "",
              "| Case | Source work [J] | Response residual / work [%] | Final peak speed [m/s] | Final kinetic energy [J] | Numerical / physical heat [%] | Additional adaptive limited stages |",
              "|---|---:|---:|---:|---:|---:|---:|"]
    for r in m["pulse_summary"]:
        d, s = r["diagnostics"], r["samples"][-1]
        lines.append(f"| {r['name']} | {d['source_work_j']:.9g} | {100*d['response_residual_over_source_work']:.7g} | {s['maximum_speed_m_s']:.7g} | {s['kinetic_energy_j']:.7g} | {100*d['numerical_interface_heat_j']/d['ohmic_heat_j']:.6g} | {d['limited_euler_stages']} |")
    lines += ["", "Source work excludes external-field self/mutual and coil-circuit energy. Numerical interface heat is separate from physical Ohmic heat and does not include all implicit Rusanov/UCT dissipation. The characteristic rise diffusion length √(D t_rise)=0.0632456 m spans about 1.01, 2.02, 4.05 and 8.10 radial cells on these grids; this is a scale estimate, not the exact finite-coil skin profile.", "",
              "| Method at nr=128 | Time [s] | Density min / max [kg/m³] | J²-weighted RMS depth [m] | Depth / radial cell | Max local numerical / physical eta [%] |",
              "|---|---:|---:|---:|---:|---:|"]
    for r in m["pulse_summary"]:
        if r["name"] not in ("standard_128", "strip_128"):
            continue
        for s in r["samples"]:
            low, high = s["density_min_max_kg_m3"]
            lines.append(f"| {r['method']} | {s['time_s']:.3g} | {low:.6g} / {high:.6g} | {s['current_weighted_rms_depth_m']:.6g} | {128*s['current_weighted_rms_depth_m']:.5g} | {100*s['maximum_numerical_eta_over_physical']:.6g} |")
    lines += ["", "Current depth is measured from the nearest gas wall/end with cylindrical nodal J² weights. It is a width diagnostic rather than an assumed exponential fit. Local eta ratios and integrated heat ratios measure different aspects of numerical dissipation.", "",
              "### Full-field spatial differences", "",
              "Cell fields are restricted from the finer mesh by exact cylindrical volume weights; currents use coincident nodes with cylindrical dual weights. Values below are the largest relative L2 difference across the three saved times. Density/pressure/B denominators use perturbations from the initial gas/bias; velocity/current use their reference norms. These are adjacent-grid differences, not analytic errors or proof of the finest grid's absolute accuracy.", "",
              "| Coarse / fine | Velocity [%] | Density perturbation [%] | Pressure perturbation [%] | Magnetic perturbation [%] | Current [%] |",
              "|---|---:|---:|---:|---:|---:|"]
    for r in m["pulse_refinement"]:
        values = [100*max(t[key] for t in r["history"]) for key in
                  ("velocity_relative_L2", "density_perturbation_relative_L2", "pressure_perturbation_relative_L2", "magnetic_perturbation_relative_L2", "current_relative_L2")]
        lines.append(f"| {r['coarse_or_test']} / {r['reference']} | "+" | ".join(f"{v:.6g}" for v in values)+" |")
    lines += ["", "### Same-grid method, time and matching controls", "",
              "The time control halves both CFL and the prescribed-source time cap. The farther matching control preserves the same actual gas domain and sources.", "",
              "| Test / reference | Max velocity L2 difference [%] | Max magnetic perturbation L2 difference [%] | Max current L2 difference [%] |",
              "|---|---:|---:|---:|"]
    for r in m["pulse_controls"]:
        values = [100*max(t[key] for t in r["history"]) for key in
                  ("velocity_relative_L2", "magnetic_perturbation_relative_L2", "current_relative_L2")]
        lines.append(f"| {r['coarse_or_test']} / {r['reference']} | "+" | ".join(f"{v:.6g}" for v in values)+" |")
    lines += ["", "## Gas-end sensitivity with a common magnetic boundary", "",
              "At nr=64 and matching extent 3 m, move reflecting gas ends from z=±1 to ±2 m with identical coils, initial gas state and spacing. The physical gas domain and total mass change. Compare common cells at |z|<0.5 m; source work refers to the whole actual domain. This quantifies a boundary-condition choice rather than demonstrating an open exhaust or free surface.", "",
              "| Method | Short/long source-work difference [% of long] | Max central velocity difference [%] | Max central magnetic perturbation difference [%] |",
              "|---|---:|---:|---:|"]
    for r in m["end_controls"]:
        lines.append(f"| {r['short_case'].split('_')[0]} | {100*r['source_work_change_relative_to_long']:.6g} | {100*max(t['velocity_relative_L2'] for t in r['history']):.6g} | {100*max(t['magnetic_perturbation_relative_L2'] for t in r['history']):.6g} |")
    method_fine = next(r for r in m["pulse_controls"] if r["coarse_or_test"] == "strip_128")
    flow_difference = 100*max(t["velocity_relative_L2"] for t in method_fine["history"])
    standard_fine = next(r for r in m["pulse_refinement"] if r["coarse_or_test"] == "standard_64")
    strip_fine = next(r for r in m["pulse_refinement"] if r["coarse_or_test"] == "strip_64")
    standard_change = 100*max(t["velocity_relative_L2"] for t in standard_fine["history"])
    strip_change = 100*max(t["velocity_relative_L2"] for t in strip_fine["history"])
    lines += ["", "## Interpretation and remaining scope", "",
              f"The strip has a material accuracy cost in the driven gas: at nr=128 the maximum same-grid velocity L2 difference is {flow_difference:.4g}%. The 64-to-128 changes are {standard_change:.4g}% for the standard method and {strip_change:.4g}% for the strip. These measured differences prevent treating the strip's small shock energy residual as evidence of an accurate general-purpose boundary replacement. It stays experimental and no production default is promoted.", "",
              "The fixed strip is a controlled numerical option, not material physics. Assess its shock energy and analytic skin results together with driven-flow grid differences and the changed-end controls. A positive state, small response residual, or close agreement between two methods cannot establish a physically accurate formation experiment.", "",
              "The evolving free-plasma/vacuum interface, actual gas loss/ends, regime-appropriate transport and ionization, broader shock/cold-flow and curved-equilibrium accuracy, converged formation and dynamic experiment UI remain open. No exterior fluid atmosphere, energy correction or transport closure is introduced by this study. Stage 2 formation remains gated.", "",
              "Full inputs, independent references, profiles, field hashes and source/environment records are in `outputs/boundary_verification/metrics.json` and its case directories. The original compatible-energy records retain their own provenance. Rebuild with `python -m examples.verify_boundary` after all cases finish.", "",
              "![Boundary accuracy study](../outputs/boundary_verification/boundary_verification.png)"]
    Path("docs/BOUNDARY_VERIFICATION.md").write_text("\n".join(lines)+"\n")


def verify(directory="outputs/boundary_verification"):
    directory = Path(directory); source = source_record()
    pulses = {name: read_case(directory/"pulses", name) for name, *_ in protocol()}
    if any(row[0]["provenance"] != source for row in pulses.values()):
        raise ValueError("Pulse source does not match the current frozen study")
    skin = []
    for method in ("physical", "standard", "strip"):
        for n, cfl in [(n, .2) for n in (16, 32, 64, 128)]+[(64, .1)]:
            row = json.loads((directory/"skin"/f"{method}_{n}_cfl{cfl}.json").read_text())
            for name in ("cylindrical_skin_reference", "cylindrical_skin_measurements"):
                if row["provenance"][name+"_sha256"] != hashlib.sha256(Path(f"examples/{name}.py").read_bytes()).hexdigest():
                    raise ValueError("Skin case source mismatch")
            skin.append(row)
    shocks = {"standard": [], "strip": []}
    old = Path("outputs/compatible_energy_verification")
    for n in (24, 48, 96, 192, 384):
        shocks["standard"].append(json.loads((old/"cases"/f"shock_{n}.json").read_text()))
        path = (old/"boundary_order_probe" if n < 384 else directory/"fine_shock")/f"shock_{n}_cfl0.2.json"
        shocks["strip"].append(json.loads(path.read_text()))
    ends = {f"{method}_{length}": read_case(directory/"ends", f"{method}_length{length}_64")
            for method in ("standard", "strip") for length in (1, 2)}
    for row, _ in ends.values():
        expected = {**source, "end_controls_sha256": hashlib.sha256(Path("examples/boundary_end_controls.py").read_bytes()).hexdigest()}
        if row["provenance"] != expected:
            raise ValueError("End-control source mismatch")
    all_rows = skin+[r for rows in shocks.values() for r in rows]+[r[0] for r in ends.values()]
    common_keys = ("source_sha256", "python", "platform", "versions", "interpretation")
    if any(not r["completed"] or any(r["provenance"][key] != source[key] for key in common_keys) for r in all_rows):
        raise ValueError("Incomplete or mismatched-source study record")
    if any(r["provenance"]["boundary_probe_sha256"] != source["boundary_order_probe_sha256"] for r in shocks["strip"]):
        raise ValueError("Shock strip source mismatch")
    m = {"provenance": source, "shocks": shocks, "skin": skin,
         "pulse_records": [r[0] for r in pulses.values()], "end_records": [r[0] for r in ends.values()],
         "pulse_summary": [summarize(case) for case in pulses.values()],
         "pulse_refinement": [compare(pulses[f"{method}_{n}"], pulses[f"{method}_{2*n}"])
                              for method in ("standard", "strip") for n in (16, 32, 64)],
         "pulse_controls": [compare(pulses[f"strip_{n}"], pulses[f"standard_{n}"]) for n in (16, 32, 64, 128)],
         "end_controls": [compare_ends(ends[f"{method}_1"], ends[f"{method}_2"]) for method in ("standard", "strip")]}
    for method in ("standard", "strip"):
        m["pulse_controls"] += [compare(pulses[f"{method}_64"], pulses[f"{method}_half_time"]),
                                compare(pulses[f"{method}_32"], pulses[f"{method}_far_matching"]),
                                compare(pulses[f"{method}_64"], ends[f"{method}_1"])]
    m["report_sources_sha256"] = {name: hashlib.sha256(Path(f"examples/{name}.py").read_bytes()).hexdigest()
                                  for name in ("verify_boundary", "boundary_analysis", "boundary_plot", "verification_records")}
    write_record(directory/"metrics.json", m); report(m); plot(m, directory)
    return m


if __name__ == "__main__":
    verify()
