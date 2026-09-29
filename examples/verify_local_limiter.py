"""Compare local/global accuracy, limiter support and independent energy errors."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from .local_limiter_measurements import run
from .audit_interface_energy import run as run_energy_audit
from .verification_records import write_record


def comparisons(metrics):
    rows = {r["name"]: r for r in metrics["cases"]}
    result = {}
    for name, first, second in (("shock_time", "shock_half_cfl", "shock_local_192"),
                                ("shock_ends", "shock_long_chamber", "shock_local_96")):
        a, b = rows[first], rows[second]
        entry = {"both_completed": a["completed"] and b["completed"]}
        if entry["both_completed"]:
            select = (np.asarray(a["profile"]["r_m"]) > .44)&(np.asarray(a["profile"]["r_m"]) < .66)
            entry["relative_profile_L2"] = {
                key: float(np.linalg.norm(np.asarray(a["profile"][key])[select]-np.asarray(b["profile"][key])[select])/
                           np.linalg.norm(np.asarray(b["profile"][key])[select]))
                for key in ("density", "pressure", "Bz_over_upstream_B")}
        result[name] = entry
    for name, first, second in (("pulse_time", "pulse_half_cfl", "pulse_32"),
                                ("pulse_matching", "pulse_far_matching", "pulse_16")):
        a, b = rows[first], rows[second]
        result[name] = {"both_completed": a["completed"] and b["completed"],
                        "relative_source_work_change": a["diagnostics"]["source_work_j"]/b["diagnostics"]["source_work_j"]-1,
                        "relative_peak_speed_change": a["diagnostics"]["maximum_speed_m_s"]/b["diagnostics"]["maximum_speed_m_s"]-1}
    historical = json.loads(Path("outputs/conservative_verification/metrics.json").read_text())
    old = {r["name"]: r for r in historical["cases"]}
    result["historical_global_source"] = historical["provenance"]
    result["inactive_limiter_regression"] = {
        name: {"completed": r["completed"], "limited_stages": r["diagnostics"]["limited_euler_stages"],
               "response_energy_change_j": r["diagnostics"]["response_energy_j"]-old[name]["diagnostics"]["response_energy_j"],
               "source_work_change_j": r["diagnostics"]["source_work_j"]-old[name]["diagnostics"]["source_work_j"],
               "smooth_error_change": r.get("flux_relative_L2", 0)-old[name].get("flux_relative_L2", 0)}
        for name, r in rows.items() if r["kind"] != "shock"}
    return result


def plot(metrics, directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    from plasma_bench.visualization import STYLE
    rows = {r["name"]: r for r in metrics["cases"]}
    colors = {"local": "#72ddc4", "global": "#ffcc75"}
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(16, 9.5), layout="constrained")
        for mode, cells in (("local", (24, 48, 96, 192, 384)), ("global", (96, 192))):
            selected = [rows[f"shock_{mode}_{n}"] for n in cells
                        if rows.get(f"shock_{mode}_{n}", {}).get("completed", False)]
            x = [r["radial_cells"] for r in selected]
            axes[0, 0].loglog(x, [r["errors"]["Bz_L1_over_upstream_B"] for r in selected],
                              "o-", color=colors[mode], label=f"{mode.capitalize()} limiting")
            axes[1, 1].plot(x, [r["diagnostics"]["response_energy_balance_relative"] for r in selected],
                            "o-", color=colors[mode], label=f"{mode.capitalize()} limiting")
            p = rows[f"shock_{mode}_192"]["profile"]
            axes[0, 1].plot(p["r_m"], p["Bz_over_upstream_B"], color=colors[mode], label=f"{mode.capitalize()} · nr=192")
            axes[0, 2].plot(p["r_m"], p["pressure"], color=colors[mode], label=f"{mode.capitalize()} · nr=192")
        p = rows["shock_local_192"]["profile"]
        axes[0, 1].plot(p["r_m"], p["exact_Bz_over_upstream_B"], "--", color="#ccd7e8", label="Independent reference")
        axes[0, 2].plot(p["r_m"], p["exact_pressure"], "--", color="#ccd7e8", label="Independent reference")
        axes[0, 0].set(title="Original shock · magnetic profile error", xlabel="Radial cells", ylabel="Bz L1 / upstream field")
        axes[0, 0].set_xticks([24,48,96,192,384], labels=["24","48","96","192","384"])
        axes[0, 0].xaxis.set_minor_formatter(NullFormatter())
        axes[0, 1].set(title="Interior accuracy restored locally", xlabel="Radius [m]", ylabel="Bz / upstream field", xlim=(.44,.66))
        axes[0, 2].set(title="Pressure still has a finite shock error", xlabel="Radius [m]", ylabel="Pressure [Pa]", xlim=(.44,.66))
        row = rows["shock_local_192"]; data = row["limiter_map"]
        fraction = np.zeros(data["shape"])
        if data["indices"]:
            i, j = np.asarray(data["indices"]).T
            fraction[i, j] = 100*np.asarray(data["affected_stage_counts"])/(2*row["diagnostics"]["steps"])
        mesh = row["fluid_mesh"]
        im = axes[1, 0].imshow(fraction, origin="lower", aspect="auto",
                               extent=(-mesh["half_length"], mesh["half_length"], 0, mesh["radius"]),
                               cmap="magma", interpolation="nearest")
        fig.colorbar(im, ax=axes[1, 0], label="Affected Euler stages per cell [%]", shrink=.85)
        axes[1, 0].plot([-.1,.1,.1,-.1,-.1], [.44,.44,.66,.66,.44], color="#72ddc4", linewidth=1)
        axes[1, 0].set(title="nr=192 · recorded limiter support", xlabel="Axial position [m]", ylabel="Radius [m]")
        axes[1, 1].axhline(0, color="#ccd7e8", linewidth=.6)
        axes[1, 1].set(title="Whole-space energy error is not forced exact", xlabel="Radial cells", ylabel="Signed response residual / initial energy")
        pulses = [rows[f"pulse_{n}"] for n in (16,32,64)]
        axes[1, 2].loglog([r["radial_cells"] for r in pulses],
                          [100*abs(r["diagnostics"]["response_residual_over_source_work"]) for r in pulses],
                          "o-", color="#aa9cfc", label="Local method · limiter inactive")
        axes[1, 2].set(title="Coil pulse · independent source budget", xlabel="Radial cells", ylabel="Absolute response error / applied work [%]")
        axes[1, 2].set_xticks([16,32,64], labels=["16","32","64"])
        axes[1, 2].xaxis.set_minor_formatter(NullFormatter())
        for ax in axes.flat:
            if ax != axes[1, 0]:
                ax.grid(alpha=.2); ax.legend(fontsize=7, labelcolor="#c9d5e5")
        fig.suptitle("LOCAL CONSERVATIVE / CT FALLBACK\nOriginal sharp shock · shared fluxes · no cell floors · separate vacuum-energy accuracy check",
                     fontsize=15, weight="bold")
        fig.savefig(directory/"local_limiter.png", dpi=150); plt.close(fig)


def report(metrics):
    rows = {r["name"]: r for r in metrics["cases"]}
    local, global_ = rows["shock_local_192"], rows["shock_global_192"]
    le, ge = local["errors"], global_["errors"]
    lines = ["# Local conservative/CT limiter verification", "",
             "Generated by `python -m examples.verify_local_limiter`. See [method, equations and boundary contracts](CONSERVATIVE_COUPLING.md). Historical global and split reports remain separate.", "",
             "The default complete-update limiter now replaces high-order face fluxes, corner electric fields and cylindrical metric sources only around inadmissible cells. Shared entities retain conservation and CT; newly troubled neighbors are rechecked. A nonphysical first-order baseline still fails without floors or energy repair. This is local binary order reduction, not an unconditional positivity theorem.", "",
             "## Original independent cylindrical shock", "",
             "The sharp initial data, independent Rankine–Hugoniot/cylindrical steady-flow solution, r=0.44–0.66 m and |z|<=0.1 m error window, physical eta=0 and requested 0.02 s final time are unchanged. Numerical interface diffusion remains explicit. No warm buffer or easier replacement boundary fixture was used.", "",
             "| Case | Finished | Density L1 | Pressure L1 [Pa] | Bz L1/B1 | Frozen B/rho defect | Min Bz/B1 |",
             "|---|---|---:|---:|---:|---:|---:|"]
    for mode, cells in (("local", (24,48,96,192,384)), ("global", (96,192))):
        for n in cells:
            r = rows[f"shock_{mode}_{n}"]; e = r["errors"]
            lines.append(f"| {mode} {n} | {r['completed']} | {e['density_L1_kg_m3']:.7g} | {e['pressure_L1_pa']:.7g} | {e['Bz_L1_over_upstream_B']:.7g} | {e['frozen_ratio_defect_L1']:.5g} | {e['minimum_Bz_over_upstream_B']:.7g} |")
    lines += ["", f"At nr=192, local replacement reduces magnetic L1 error by {1-le['Bz_L1_over_upstream_B']/ge['Bz_L1_over_upstream_B']:.3%} and pressure L1 error by {1-le['pressure_L1_pa']/ge['pressure_L1_pa']:.3%} relative to the freshly measured global control. A sharper interior solution does **not** establish overall physical validity or free-plasma formation.", "",
              "## Locality and complete-update robustness", "",
              "| Case | Limited / all Euler stages | Affected cell-stages [%] | Max cells affected in one stage [%] | Max passes | Max depth from gas boundary [m] | Directly affected reference-window cell-stages |",
              "|---|---:|---:|---:|---:|---:|---:|"]
    for n in (24,48,96,192,384):
        r = rows[f"shock_local_{n}"]; d, m = r["diagnostics"], r["limiter_map"]
        fraction = d["affected_euler_cell_stages"]/max(1, 2*d["steps"]*np.prod(m["shape"]))
        lines.append(f"| local {n} | {d['limited_euler_stages']} / {2*d['steps']} | {100*fraction:.6g} | {100*d['maximum_limiter_cell_fraction']:.6g} | {d['maximum_limiter_iterations']} | {m['maximum_depth_from_gas_boundary_m']:.6g} | {m['affected_reference_window_cell_stages']} |")
    lines += ["", "Affected cells include diagonal neighbors sharing a replaced corner. Untouched update stencils remain high order; subsequent physical/numerical propagation can still carry boundary errors inward. Minimum high-order weight is zero in a locally replaced stencil, one elsewhere. The plotted box is the independent shock error window.", "",
              "Manufactured tests independently check an energy-transfer chain that requires repeated neighbor detection, a magnetic-energy-only inadmissibility, one-dimensional Bz/rho preservation, untouched distant derivatives, boundary energy exchange, CT and invalid-baseline failure. The original nr=96 shock also exercises accepted local limiter maps in the test suite.", "",
              "## Independent energy budgets", "",
              "| Case | Fluid energy-flux residual [J] | Response residual [J] | Response residual / initial energy | Net outward fluid energy [J] | Numerical interface heat [J] |",
              "|---|---:|---:|---:|---:|---:|"]
    for mode, cells in (("local", (24,48,96,192,384)), ("global", (96,192))):
        for n in cells:
            r = rows[f"shock_{mode}_{n}"]; d = r["diagnostics"]
            lines.append(f"| {mode} {n} | {d['fluid_energy_balance_residual_j']:.5g} | {d['response_energy_balance_residual_j']:.6g} | {d['response_energy_balance_relative']:.6g} | {d['outward_fluid_energy_j']:.7g} | {d['numerical_interface_heat_j']:.7g} |")
    ld, gd = local["diagnostics"], global_["diagnostics"]
    fine = rows["shock_local_384"]; fd = fine["diagnostics"]
    lines += ["", f"At nr=192 the signed full-space response error is {ld['response_energy_balance_relative']:.6g} locally versus {gd['response_energy_balance_relative']:.6g} globally. Its absolute value increases from local nr=96 to nr=192 despite improved central shock profiles. Exact fluid-cell conservation does not remove mismatch among cell magnetic energy, nodal full-space magnetic energy, boundary traces and selected CT/energy fluxes. This unresolved boundary/energy accuracy issue is retained explicitly; no energy correction forces agreement.", "",
              f"The local nr=192 absolute residual {ld['response_energy_balance_residual_j']:.7g} J is {abs(ld['response_energy_balance_residual_j']/ld['outward_fluid_energy_j']):.4%} of the net outward fluid energy, an additional transfer scale. That boundary transfer is not external source work, and normalization by the much larger initial energy alone would understate this interface issue.", "",
              f"The nr=384 case completes with central B error {fine['errors']['Bz_L1_over_upstream_B']:.7g} and pressure error {fine['errors']['pressure_L1_pa']:.7g} Pa, but the response residual grows to {fd['response_energy_balance_residual_j']:.7g} J, or {abs(fd['response_energy_balance_residual_j']/fd['outward_fluid_energy_j']):.4%} of net outward fluid energy. Numerical interface heat also rises to {fd['numerical_interface_heat_j']:.7g} J. This refinement does **not** converge the interface energy budget. The minimum density is {fd['minimum_density_kg_m3']:.7g} kg/m³, without a density floor; the run takes {fd['steps']} steps, {fine['setup_seconds']:.1f} s setup and {fine['advance_seconds']:.1f} s advancement on the recorded CPU environment. Limiter support remains outside the central reference window. Strong wall-layer dynamics and their numerical/physical validity remain unresolved.", "",
              "The response budget is gas plus induced self-field energy through the vacuum to infinity. External self/mutual and circuit energy are excluded. Shock source work and physical Ohmic heat are zero; interface heat is numerical, and does not include all bulk Rusanov dissipation.", "",
              "## Instantaneous energy-rate audit", "",
              "The [independently derived rate identity](INTERFACE_ENERGY_AUDIT.md) separates outward fluid power, cell magnetic-energy rate and nodal magnetic work. High, low and locally selected rates below are evaluated at the same accepted state. A symmetric perturbation recompletes the magnetic vacuum and differences the integrated response energy at three sizes. The smallest absolute identity error is reported; all values are retained in the JSON. This checks the diagnostic identity, not the physical accuracy of the fluxes.", "",
              "| Cells | Time [s] | High-order defect rate [W] | Low-order rate [W] | Selected rate [W] | Selected minus high [W] | Best finite-difference identity error [W] |",
              "|---:|---:|---:|---:|---:|---:|---:|"]
    for audit in metrics["energy_audit"]:
        for s in audit["snapshots"]:
            rates = s["rates"]
            lines.append(f"| {audit['radial_cells']} | {s['time_s']:.4g} | {rates['high']['response_defect_rate_w']:.7g} | {rates['low']['response_defect_rate_w']:.7g} | {rates['selected']['response_defect_rate_w']:.7g} | {s['selected_minus_high_w']:.7g} | {min(d['absolute_identity_error_w'] for d in s['finite_differences']):.5g} |")
    lines += ["", "A nonzero defect rate exists before time integration. A small selected-minus-high difference at a snapshot does not remove effects of earlier fallback on the evolving wall state. Spatial consistency, nonlinear boundary evolution and limiter/time sensitivity must be distinguished before changing the scheme.", "",
              "## Controls and retained failure records", "",
              "| Case | Finished | Accepted time [s] | Failure if any |", "|---|---|---:|---|"]
    for r in rows.values():
        if not r["completed"] or r["name"] in ("shock_half_cfl", "shock_long_chamber"):
            lines.append(f"| {r['name']} | {r['completed']} | {r['diagnostics']['time_s']:.8g} | {r.get('failure', '—')} |")
    lines += ["", "Any failed row reports its last accepted state, not a final-time convergence point. Both time and end-location comparisons use a resolution at which local limiting is active.", "",
              "| Shock control | Both completed | Density relative profile L2 change | Pressure change | Bz change |",
              "|---|---|---:|---:|---:|"]
    for name in ("shock_time", "shock_ends"):
        c = metrics["comparisons"][name]; p = c.get("relative_profile_L2")
        values = " | ".join(f"{p[k]:.6%}" for k in ("density", "pressure", "Bz_over_upstream_B")) if p else "— | — | —"
        lines.append(f"| {name} | {c['both_completed']} | {values} |")
    lines += ["", "| Shock control | Response residual [J] | Relative response residual | Max affected cells per stage [%] |",
              "|---|---:|---:|---:|"]
    for name in ("shock_local_192", "shock_half_cfl", "shock_local_96", "shock_long_chamber"):
        d = rows[name]["diagnostics"]
        lines.append(f"| {name} | {d['response_energy_balance_residual_j']:.7g} | {d['response_energy_balance_relative']:.7g} | {100*d['maximum_limiter_cell_fraction']:.6g} |")
    lines += ["", "## Coil pulse and smooth-solution regression", "",
              "The same prescribed-coil pulse and analytic advection/diffusion fixtures are rerun with the local default. The historical comparison fingerprint is retained in the JSON. Inactive limiter regression reports changes in source work, response energy and analytic error, rather than assuming unchanged behavior.", "",
              "| Pulse case | Work [J] | Response error / work [%] | Physical / numerical interface heat [J] | Peak speed [m/s] | Limited stages |",
              "|---|---:|---:|---:|---:|---:|"]
    for name in ("pulse_16", "pulse_32", "pulse_64", "pulse_half_cfl", "pulse_far_matching"):
        d = rows[name]["diagnostics"]
        lines.append(f"| {name} | {d['source_work_j']:.8g} | {100*d['response_residual_over_source_work']:.7g} | {d['ohmic_heat_j']:.7g} / {d['numerical_interface_heat_j']:.7g} | {d['maximum_speed_m_s']:.7g} | {d['limited_euler_stages']} |")
    lines += ["", "| Pulse control | Relative work change | Relative peak-speed change |", "|---|---:|---:|"]
    for name in ("pulse_time", "pulse_matching"):
        c = metrics["comparisons"][name]
        lines.append(f"| {name} | {c['relative_source_work_change']:.6%} | {c['relative_peak_speed_change']:.6%} |")
    lines += ["", "| Smooth case | Flux relative L2 | Limited stages |", "|---|---:|---:|"]
    for kind in ("advection", "diffusion"):
        for n in (12,24,48):
            r = rows[f"{kind}_{n}"]
            lines.append(f"| {kind} {n} | {r['flux_relative_L2']:.8g} | {r['diagnostics']['limited_euler_stages']} |")
    regression = metrics["comparisons"]["inactive_limiter_regression"]
    lines += ["", f"Maximum absolute change from historical inactive-limiter cases: source work {max(abs(r['source_work_change_j']) for r in regression.values()):.5g} J, response energy {max(abs(r['response_energy_change_j']) for r in regression.values()):.5g} J, smooth relative L2 error {max(abs(r['smooth_error_change']) for r in regression.values()):.5g}. These controls retain their prior limitations: driven peak flow and interface heat still require refinement.", ""]
    completed = [r for r in rows.values() if r["completed"]]
    lines += [f"Across completed cases, maximum absolute relative mass change is {max(abs(r['diagnostics']['mass_relative_change']) for r in completed):.5g}; compatible face divergence {max(r['diagnostics']['max_face_divergence_t_m'] for r in completed):.5g} T/m; vacuum relative residual {max(r['diagnostics']['vacuum_relative_residual'] for r in completed):.5g}. Complete inputs, histories, sparse limiter maps, radial profiles, accepted-time failures and source/environment fingerprints are in `outputs/local_limiter_verification/metrics.json`.", "",
              "![Local limiter accuracy and locality](../outputs/local_limiter_verification/local_limiter.png)", "",
              "**Remaining scope:** reconcile interface/CT and full-space energy accuracy without thermal repairs, establish broader shock and cold-rarefaction accuracy, resolve driven skin layers, implement a defensible moving plasma/vacuum boundary and actual ends, add regime-appropriate transport/ionization/Hall physics, and perform converged formation experiments with a dynamic experimental interface. The full FRC goal remains active; the original Stage 2 gate stays closed."]
    Path("docs/LOCAL_LIMITER_VERIFICATION.md").write_text("\n".join(lines)+"\n")


def verify(directory="outputs/local_limiter_verification", *, reuse=False):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    metrics = run(directory/"cases", reuse=reuse)
    metrics["comparisons"] = comparisons(metrics)
    metrics["energy_audit"] = []
    audit_hash = hashlib.sha256(Path("examples/audit_interface_energy.py").read_bytes()).hexdigest()
    for n in (96,192):
        path = directory/"energy_audit"/f"shock_{n}.json"
        audit = json.loads(path.read_text()) if reuse and path.exists() else None
        if (not audit or len(audit["snapshots"]) != 5 or audit["snapshots"][-1]["time_s"] != .02 or
            audit["audit_source_sha256"] != audit_hash or
            audit["provenance"]["source_sha256"] != metrics["provenance"]["source_sha256"]):
            audit = run_energy_audit(n, directory/"energy_audit")
        metrics["energy_audit"].append(audit)
    metrics["report_source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    metrics["serialization_source_sha256"] = hashlib.sha256(Path("examples/verification_records.py").read_bytes()).hexdigest()
    for row in metrics["cases"]:
        write_record(directory/"cases"/f"{row['name']}.json", row)
    write_record(directory/"metrics.json", metrics)
    plot(metrics, directory); report(metrics)
    return metrics


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse", action="store_true")
    verify(reuse=parser.parse_args().reuse)
