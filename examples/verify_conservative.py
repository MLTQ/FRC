"""Report conservative shock accuracy, boundary limitations and independent budgets."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from .conservative_measurements import all_measurements


def comparisons(metrics):
    rows = {r["name"]: r for r in metrics["cases"]}
    result = {}
    for name, first, second in (("shock_time", "shock_half_cfl", "shock_96"),
                                ("shock_ends", "shock_long_chamber", "shock_48")):
        a, b = rows[first], rows[second]
        select = (np.asarray(a["profile"]["r_m"]) > .44)&(np.asarray(a["profile"]["r_m"]) < .66)
        values = {}
        for key in ("density", "pressure", "Bz_over_upstream_B"):
            x, y = np.asarray(a["profile"][key])[select], np.asarray(b["profile"][key])[select]
            values[key] = float(np.linalg.norm(x-y)/np.linalg.norm(y))
        result[name] = {"both_completed": a["completed"] and b["completed"], "relative_profile_L2": values}
    for name, first, second in (("pulse_time", "pulse_half_cfl", "pulse_32"),
                                ("pulse_matching", "pulse_far_matching", "pulse_16"),
                                ("pulse_interface", "pulse_no_interface_diffusion", "pulse_32")):
        a, b = rows[first], rows[second]
        result[name] = {"both_completed": a["completed"] and b["completed"],
                        "relative_source_work_change": a["diagnostics"]["source_work_j"]/b["diagnostics"]["source_work_j"]-1,
                        "relative_peak_speed_change": a["diagnostics"]["maximum_speed_m_s"]/b["diagnostics"]["maximum_speed_m_s"]-1}
    return result


def plot(metrics, directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    from plasma_bench.visualization import STYLE
    rows = {r["name"]: r for r in metrics["cases"]}
    shocks = [rows[f"shock_{n}"] for n in (24, 48, 96, 192) if rows[f"shock_{n}"]["completed"]]
    pulses = [rows[f"pulse_{n}"] for n in (8, 16, 32, 64) if rows[f"pulse_{n}"]["completed"]]
    historical = json.loads(Path("outputs/shock_verification/metrics.json").read_text())
    old = {r["name"]: r for r in historical["cases"]}
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(16, 9), layout="constrained")
        for row, label, color in ((old["centered_48"], "Earlier centered split", "#f687a6"),
                                  (old["diffusive_48"], "Earlier split + volume diffusion", "#ffcc75"),
                                  (rows["shock_48"], "Conservative + magnetic vacuum", "#72ddc4")):
            p = row["profile"]
            axes[0, 0].plot(p["r_m"], np.asarray(p["Bz_over_upstream_B"])-p["density"], ".-", color=color, label=label)
        axes[0, 0].set(xlim=(.44, .66), title="48 cells · frozen-flux defect", xlabel="Radius [m]", ylabel="Bz/B1 − rho/rho1")
        for selected, label, color in ((shocks, "Conservative + vacuum", "#72ddc4"),
                                      ([old[f"diffusive_{n}"] for n in (24,48,96,192)], "Earlier split + volume diffusion", "#ffcc75")):
            axes[0, 1].loglog([r["radial_cells"] for r in selected], [r["errors"]["Bz_L1_over_upstream_B"] for r in selected], "o-", color=color, label=label)
        axes[0, 1].set(title="Independent shock · magnetic profile error", xlabel="Radial cells", ylabel="Bz L1 error / upstream B")
        axes[0, 2].loglog([r["radial_cells"] for r in shocks], [abs(r["diagnostics"]["response_energy_balance_relative"]) for r in shocks], "o-", color="#aa9cfc", label="Full-space response error / initial energy")
        axes[0, 2].set(title="Energy error measured, not corrected", xlabel="Radial cells", ylabel="Absolute relative response error")
        axes[1, 0].loglog([r["radial_cells"] for r in pulses], [100*abs(r["diagnostics"]["response_residual_over_source_work"]) for r in pulses], "o-", color="#ffcc75", label="Final error relative to actual coil work")
        axes[1, 0].set(title="Driven pulse · work-budget refinement", xlabel="Radial cells", ylabel="Response error / source work [%]")
        for kind, color in (("advection", "#72ddc4"), ("diffusion", "#aa9cfc")):
            selected = [rows[f"{kind}_{n}"] for n in (12,24,48)]
            axes[1, 1].loglog([r["radial_cells"] for r in selected], [r["flux_relative_L2"] for r in selected], "o-", color=color, label=kind.capitalize())
        axes[1, 1].set(title="Independent smooth-solution errors", xlabel="Radial cells", ylabel="Flux relative L2")
        axes[1, 2].plot([r["radial_cells"] for r in pulses], [r["diagnostics"]["maximum_speed_m_s"] for r in pulses], "o-", color="#72ddc4", label="Conservative pulse · final peak speed")
        axes[1, 2].set(title="Driven-flow convergence remains a requirement", xlabel="Radial cells", ylabel="Peak speed [m/s]")
        for ax in axes.flat:
            ax.grid(alpha=.2); ax.legend(fontsize=7, labelcolor="#c9d5e5")
        for ax in axes[0, 1:]:
            ax.set_xticks([24,48,96,192], labels=["24","48","96","192"]); ax.xaxis.set_minor_formatter(NullFormatter())
        axes[1, 0].set_xticks([8,16,32,64], labels=["8","16","32","64"]); axes[1, 0].xaxis.set_minor_formatter(NullFormatter())
        axes[1, 1].set_xticks([12,24,48], labels=["12","24","48"]); axes[1, 1].xaxis.set_minor_formatter(NullFormatter())
        fig.suptitle("CONSERVATIVE PLASMA + MAGNETIC VACUUM\nShared shock fluxes · explicit interface diffusion · complete-update limiter activity recorded", fontsize=15, weight="bold")
        fig.savefig(directory/"conservative.png", dpi=150); plt.close(fig)


def report(metrics):
    rows = {r["name"]: r for r in metrics["cases"]}
    finest = rows["shock_192"]
    old_rows = {r["name"]: r for r in json.loads(Path("outputs/shock_verification/metrics.json").read_text())["cases"]}
    old_finest = old_rows["diffusive_192"]
    lines = ["# Conservative plasma/vacuum verification", "",
             "Generated by `python -m examples.verify_conservative`. See [equations, numerical diffusion, boundaries and energy contracts](CONSERVATIVE_COUPLING.md).", "",
             "The new solver advances mass, momentum, total MHD energy and constrained magnetic induction together, using consistent interior Rusanov/UCT fluxes. It retains a genuine elliptic magnetic exterior. A declared positive interface diffusion and complete conservative-update limiter handle the original shock's difficult cold boundary cells. Their heat/activity are reported; they are not physical transport. FRC formation is still unvalidated.", "",
             "## Original cylindrical magnetic shock", "",
             "The independent Rankine–Hugoniot and cylindrical steady-flow fixture, initial fields, r=0.44–0.66 m / |z|<=0.1 m error window, and requested 0.02 s end time are unchanged from the prior shock study. No warm pressure buffer, smoother replacement shock or easier remote initial state is substituted. The chamber remains reflecting for gas and transparent to the infinite magnetic exterior.", "",
             "| Cells | Finished | Density L1 | Pressure L1 [Pa] | Bz L1/B1 | Frozen B/rho defect | Min Bz/B1 | Limited Euler stages | Min high-order weight |", "|---:|---|---:|---:|---:|---:|---:|---:|---:|"]
    for n in (24,48,96,192):
        r = rows[f"shock_{n}"]; e, d = r["errors"], r["diagnostics"]
        lines.append(f"| {n} | {r['completed']} | {e['density_L1_kg_m3']:.6g} | {e['pressure_L1_pa']:.6g} | {e['Bz_L1_over_upstream_B']:.6g} | {e['frozen_ratio_defect_L1']:.5g} | {e['minimum_Bz_over_upstream_B']:.6g} | {d['limited_euler_stages']} / {2*d['steps']} | {d['minimum_high_order_weight']:.5g} |")
    lines += ["", "Both mass and induction use the same interior numerical flux, removing the former split method's mismatched frozen-flux diffusion. Any remaining small defect includes actual two-dimensional influence from the chamber ends. The plotted earlier split results are historical measurements retained with their own fingerprint in `outputs/shock_verification/metrics.json`, not silently rerun or overwritten.", "",
              "The limiter uses a single global convex weight when a full Euler proposal would be nonphysical. It preserves conservative flux/CT accounting and introduces no cell floor or energy repair, but can reduce accuracy throughout the grid because of one boundary cell. It is not an assertion of uniformly second-order evolution or a theorem that every low-order plasma/vacuum proposal is admissible.", "",
              f"At nr=192, the new B-profile error is {finest['errors']['Bz_L1_over_upstream_B']:.6g} versus the earlier split C=1 error {old_finest['errors']['Bz_L1_over_upstream_B']:.6g}, and pressure error is {finest['errors']['pressure_L1_pa']:.6g} versus {old_finest['errors']['pressure_L1_pa']:.6g} Pa. The improved frozen-flux invariant is **not uniform improvement in every observable**. Localizing the flux limiter remains a numerical priority; {finest['diagnostics']['limited_euler_stages']} of {2*finest['diagnostics']['steps']} Euler stages add diffusion throughout that fine grid because of troubled boundary cells.", "",
              "| Control | Finished | Accepted time [s] | Limited stages | Failure if any |", "|---|---|---:|---:|---|"]
    for name in ("shock_no_limiter", "shock_no_interface_diffusion", "shock_first_order", "shock_half_cfl", "shock_long_chamber"):
        r = rows[name]; d = r["diagnostics"]
        lines.append(f"| {name} | {r['completed']} | {d['time_s']:.7g} | {d['limited_euler_stages']} | {r.get('failure','—')} |")
    lines += ["", "Failed controls are retained at their last accepted state and excluded from final-time convergence curves. Disabling numerical interface diffusion in its control also disables full-update limiting so the original failure remains visible.", "",
              "| Completed shock | Fluid energy-flux residual [J] | Response residual / initial energy | Numerical interface heat [J] |", "|---:|---:|---:|---:|"]
    for n in (24,48,96,192):
        r = rows[f"shock_{n}"]; d = r["diagnostics"]
        if r["completed"]:
            lines.append(f"| {n} | {d['fluid_energy_balance_residual_j']:.5g} | {d['response_energy_balance_relative']:.5g} | {d['numerical_interface_heat_j']:.6g} |")
    lines += ["", "Fluid energy flux balances to roundoff. The independent gas-plus-full-space-induced-field budget has finite discretization error and is **not repaired**. Numerical interface heating is diagnosed independently; physical eta and physical Ohmic heat are zero in this shock. Bulk Rusanov dissipation is additional numerical dissipation, not included in that interface-only heat diagnostic.", "",
              "## Prescribed-coil pulse", "",
              "The prior two-coil unit-scale SI fixture is retained: rho=1 kg/m³, p=1 Pa, rest, bias=sqrt(mu0) T, physical eta=0.1 mu0 ohm m; radius-1.23 m loops at z=±0.37 m; vacuum center amplitude=0.3 bias; rise/hold/fall=0.04/0.02/0.04 s; 0.002 s source cap. Histories have 21 output times. These are numerical verification parameters, not an argon or hardware prescription.", "",
              "| Cells / control | Work [J] | Physical Ohmic heat [J] | Numerical interface heat [J] | Final response error / work | Peak speed [m/s] | Limited stages |", "|---|---:|---:|---:|---:|---:|---:|"]
    for name in ("pulse_8", "pulse_16", "pulse_32", "pulse_64", "pulse_half_cfl", "pulse_far_matching", "pulse_no_interface_diffusion"):
        r = rows[name]; d = r["diagnostics"]
        lines.append(f"| {name} | {d['source_work_j']:.7g} | {d['ohmic_heat_j']:.7g} | {d['numerical_interface_heat_j']:.7g} | {d['response_residual_over_source_work']:.5%} | {d['maximum_speed_m_s']:.6g} | {d['limited_euler_stages']} |")
    lines += ["", "Whole-space response energy is gas energy plus induced self-field energy, including the magnetic exterior to infinity. It excludes external self/mutual and circuit energy. Relative error is normalized to actual source work here, so resting-gas energy cannot hide injection error. Numerical heat is not relabeled as physical conductivity; interface sensitivity and physical heating still require refinement.", "",
              "## Time, boundary and smooth-solution checks", "",
              "| Shock control | Density profile relative L2 change | Pressure change | Bz change |", "|---|---:|---:|---:|"]
    for name in ("shock_time", "shock_ends"):
        values = metrics["comparisons"][name]["relative_profile_L2"]
        lines.append(f"| {name} | {values['density']:.5%} | {values['pressure']:.5%} | {values['Bz_over_upstream_B']:.5%} |")
    lines += ["", "| Pulse control | Relative source-work change | Relative peak-speed change |", "|---|---:|---:|"]
    for name in ("pulse_time", "pulse_matching", "pulse_interface"):
        values = metrics["comparisons"][name]
        lines.append(f"| {name} | {values['relative_source_work_change']:.6%} | {values['relative_peak_speed_change']:.6%} |")
    lines += ["", "The long-chamber shock control doubles its half-length at unchanged dr=dz. Pulse matching moves the magnetic matching rectangle from extent 2 to 4 at nr=16. These measure specified controls, not a general free-interface guarantee. Time halving and limiter switching must be assessed together.", "",
              "| Smooth case | Cells | Flux relative L2 | Limited Euler stages |", "|---|---:|---:|---:|"]
    for kind in ("advection", "diffusion"):
        for n in (12,24,48):
            r = rows[f"{kind}_{n}"]
            lines.append(f"| {kind} | {n} | {r['flux_relative_L2']:.7g} | {r['diagnostics']['limited_euler_stages']} |")
    completed = [r for r in rows.values() if r["completed"]]
    lines += ["", f"Across completed runs, maximum relative mass change is {max(abs(r['diagnostics']['mass_relative_change']) for r in completed):.4g}; maximum compatible face divergence is {max(r['diagnostics']['max_face_divergence_t_m'] for r in completed):.4g} T/m. All case inputs, source waveforms, parameter controls, histories, profiles, failure records, versions and source fingerprint are in `outputs/conservative_verification/metrics.json`.", "",
              "![Conservative shock and vacuum-work verification](../outputs/conservative_verification/conservative.png)", "",
              "**Still required:** improve/localize flux limiting and verify stronger multidimensional shocks, validate interface traction/work and physical versus numerical dissipation for the intended regime, resolve driven skin layers, replace fixed gas boundaries with a defensible evolving plasma/vacuum treatment and actual end conditions, add transport/ionization/Hall validity, and perform converged formation studies. An initialized shock and a driven verification pulse do not demonstrate a torus, confinement or three-dimensional stability."]
    Path("docs/CONSERVATIVE_VERIFICATION.md").write_text("\n".join(lines)+"\n")


def verify(directory="outputs/conservative_verification", *, reuse=False):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    metrics = all_measurements(directory/"cases", reuse=reuse)
    metrics["comparisons"] = comparisons(metrics)
    metrics["report_source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (directory/"metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=False)+"\n")
    plot(metrics, directory); report(metrics)
    return metrics


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse", action="store_true")
    verify(reuse=parser.parse_args().reuse)
