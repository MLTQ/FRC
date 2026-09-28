"""Render and report measured shock errors without hiding centered failures."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from .shock_measurements import all_measurements


def by_name(metrics):
    return {row["name"]: row for row in metrics["cases"]}


def profile_difference(first, second):
    r = np.asarray(first["profile"]["r_m"])
    select = (r > .44)&(r < .66)
    output = {}
    for key in ("density", "pressure", "Bz_over_upstream_B"):
        a, b = np.asarray(first["profile"][key]), np.asarray(second["profile"][key])
        output[key] = float(np.linalg.norm((a-b)[select])/np.linalg.norm(b[select]))
    return output


def comparisons(metrics):
    rows = by_name(metrics)
    result = {}
    for label, candidate, baseline in (
        ("inner_extension", "inner_control", "diffusive_48"),
        ("axial_ends", "ends_control", "diffusive_48"),
        ("wall_taper_short", "wall_untapered_short", "wall_tapered_short"),
        ("CFL_0.1", "time_0.1", "diffusive_48"),
        ("CFL_0.05", "time_0.05", "time_0.1"),
        ("direct_vs_picard", "direct_comparison", "diffusive_24")):
        a, b = rows[candidate], rows[baseline]
        result[label] = {"candidate": candidate, "baseline": baseline,
                         "both_completed": a["completed"] and b["completed"],
                         "relative_profile_L2": profile_difference(a, b)}
    result["direct_vs_picard"]["advance_runtime_ratio_direct_over_picard"] = rows["direct_comparison"]["advance_seconds"]/rows["diffusive_24"]["advance_seconds"]
    return result


def plot(metrics, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    from plasma_bench.visualization import STYLE
    rows = by_name(metrics)
    styles = (("bounded", "Conservative bounded reference", "#72ddc4"),
              ("centered", "Centered magnetic split · C=0", "#f687a6"),
              ("diffusive", "Numerical diffusion · C=1", "#ffcc75"))
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(16, 9), layout="constrained")
        profile = rows["diffusive_192"]["profile"]
        axes[0, 0].plot(profile["r_m"], profile["exact_Bz_over_upstream_B"], color="white", lw=2, label="Independent steady reference")
        for prefix, label, color in styles:
            p = rows[f"{prefix}_48"]["profile"]
            axes[0, 0].plot(p["r_m"], p["Bz_over_upstream_B"], ".-", color=color, label=label)
        axes[0, 0].axhline(0, color="#b9c6d9", lw=.7, ls=":")
        axes[0, 0].set(xlim=(.44, .66), title="48 radial cells · false reversal at C=0", xlabel="Radius [m]", ylabel="Bz / upstream B")
        axes[0, 1].plot(profile["r_m"], profile["exact_pressure"], color="white", lw=2, label="Independent reference")
        for n, color in ((48, "#aa9cfc"), (96, "#72ddc4"), (192, "#ffcc75")):
            p = rows[f"diffusive_{n}"]["profile"]
            axes[0, 1].plot(p["r_m"], p["pressure"], ".-", color=color, label=f"C=1, nr={n}")
        axes[0, 1].set(xlim=(.44, .66), title="Shock heating · pressure under refinement", xlabel="Radius [m]", ylabel="Gas pressure [Pa]")
        for ax, key, title, ylabel in (
            (axes[0, 2], "Bz_L1_over_upstream_B", "Magnetic profile error", "Bz L1 error / upstream B"),
            (axes[1, 0], "pressure_L1_pa", "Pressure error", "Pressure L1 error [Pa]"),
            (axes[1, 1], "frozen_ratio_defect_L1", "Frozen B/rho defect remains", "Mean |Bz/B1 − rho/rho1|")):
            for prefix, label, color in styles:
                selected = [row for name, row in rows.items() if name.startswith(prefix+"_")]
                if prefix == "bounded" and key == "frozen_ratio_defect_L1":
                    continue  # Its roundoff values are reported in the table, not on this scale.
                ax.loglog([row["radial_cells"] for row in selected], [row["errors"][key] for row in selected], "o-", color=color, label=label)
            ax.set(title=title, xlabel="Radial cells", ylabel=ylabel)
            ax.set_xticks([24, 48, 96, 192], labels=["24", "48", "96", "192"])
            ax.xaxis.set_minor_formatter(NullFormatter())
        selected = [rows[f"diffusive_{n}"] for n in (24, 48, 96, 192)]
        axes[1, 2].plot([row["radial_cells"] for row in selected], [row["diagnostics"]["numerical_magnetic_heat_j"] for row in selected], "o-", color="#ffcc75", label="C=1 · whole chamber")
        axes[1, 2].set(title="Declared numerical magnetic heat", xlabel="Radial cells", ylabel="Cumulative numerical heat [J]")
        for ax in axes.flat:
            ax.grid(alpha=.2); ax.legend(fontsize=7, labelcolor="#c9d5e5")
        fig.suptitle("STRONG MAGNETIC SHOCK · INDEPENDENT CYLINDRICAL REFERENCE\nDiffusion reduces false magnetic structure; convergence and transport remain separate requirements", fontsize=15, weight="bold")
        fig.savefig(path, dpi=150); plt.close(fig)


def report(metrics):
    rows = by_name(metrics)
    lines = ["# Magnetic shocks: measured accuracy and remaining limitations", "",
             "Generated by `python -m examples.verify_shocks`. Equations, independent reference, boundaries and numerical-heating semantics are in [SHOCK_NUMERICS.md](SHOCK_NUMERICS.md).", "",
             "**The centered coupled method can create false magnetic reversal despite a conserved energy budget.** The explicit compressive numerical diffusion reduces this defect in the tested cases. It remains a selectable numerical regularization, off by default; it is not a physical transport model or proof of arbitrary-shock convergence. The original FRC formation entry remains gated.", "",
             "## Independent strong perpendicular shock", "",
             "At r=0.5 m, upstream rho=1 kg/m³, p=0.1 Pa, vr=-6 m/s and Bz=sqrt(mu0) T, with gamma=5/3. Independent Rankine–Hugoniot jumps give downstream rho=3.532340636 kg/m³, vr=-1.698590430 m/s, p=20.16974224 Pa and Bz/B1=3.532340636. Both cylindrical smooth branches solve their own Bernoulli invariant. The physical field remains positive. No numerical profile constructs this reference.", "",
             "Both actual two-dimensional kernels evolve independently volume-integrated initial conditions to 0.02 s. Errors below use cylindrical overlap-volume means on r=0.44–0.66 m, z=±0.1 m. The coupled fluid chamber has reflecting walls/ends and an infinite magnetic exterior; the conservative bounded comparison has periodic ends and a conducting radial wall. Remote initial extensions and boundaries are tested separately below. Pressure means the primitive observable recovered from independent conserved-variable averages, not a point sample.", "",
             "| Method | Radial cells | Density L1 [kg/m³] | Pressure L1 [Pa] | Bz L1 / B1 | Frozen B/rho defect | Minimum Bz/B1 |", "|---|---:|---:|---:|---:|---:|---:|"]
    for prefix, label in (("bounded", "Bounded conservative"), ("centered", "Centered C=0"), ("diffusive", "Numerical C=1")):
        for name, row in rows.items():
            if name.startswith(prefix+"_"):
                e = row["errors"]
                lines.append(f"| {label} | {row['radial_cells']} | {e['density_L1_kg_m3']:.6g} | {e['pressure_L1_pa']:.6g} | {e['Bz_L1_over_upstream_B']:.6g} | {e['frozen_ratio_defect_L1']:.6g} | {e['minimum_Bz_over_upstream_B']:.6g} |")
    base, corrected = rows["centered_48"]["errors"], rows["diffusive_48"]["errors"]
    lines += ["", f"At 48 radial cells, C=0 reaches Bz/B1={base['minimum_Bz_over_upstream_B']:.6g}, an unphysical sign reversal. C=1 has minimum {corrected['minimum_Bz_over_upstream_B']:.6g} and reduces magnetic L1 error by a factor {base['Bz_L1_over_upstream_B']/corrected['Bz_L1_over_upstream_B']:.3g}. This is not a general monotonicity guarantee: finite-grid overshoot and B/rho defects remain.", "",
              "The bounded conservative method preserves the special initial B/rho invariant to roundoff. The split method does not. Its C=1 frozen-ratio defect is **not monotonic under mesh refinement**, even while bulk density, pressure and B errors decrease. Gas numerical flux diffusion and the added magnetic diffusion are different. More resolution, additional shock strengths/orientations and a conservative weak-solution assessment remain necessary. A good global budget does not settle these questions.", "",
              "## Coefficient, time and boundary controls", "",
              "| Coefficient at nr=48 | Density L1 | Pressure L1 | Bz L1/B1 | Numerical magnetic heat [J] |", "|---:|---:|---:|---:|---:|"]
    for name in ("centered_48", "coefficient_0.5", "diffusive_48", "coefficient_2.0"):
        row = rows[name]; e = row["errors"]
        lines.append(f"| {row['parameters']['shock_diffusion']:g} | {e['density_L1_kg_m3']:.6g} | {e['pressure_L1_pa']:.6g} | {e['Bz_L1_over_upstream_B']:.6g} | {row['diagnostics']['numerical_magnetic_heat_j']:.6g} |")
    lines += ["", "The coefficient is not optimized into a transport claim. Reported numerical heat covers the entire chamber, including remote initial extensions. Physical resistivity and physical Ohmic heat are zero. Integrated shock dissipation need not vanish even when a local regularization coefficient vanishes with grid spacing.", "",
              "| Control | Density profile relative L2 change | Pressure change | Bz change | Both runs finished |", "|---|---:|---:|---:|---|"]
    for name, row in metrics["comparisons"].items():
        e = row["relative_profile_L2"]
        lines.append(f"| {name} | {e['density']:.5g} | {e['pressure']:.5g} | {e['Bz_over_upstream_B']:.5g} | {row['both_completed']} |")
    lines += ["", "At nr=48, inner control moves the exact-branch cutoff from r=0.38 to 0.37 m. Axial control doubles chamber half-length at fixed dr=dz. Wall control removes the velocity taper beyond r=0.85 m and compares at **0.01 s**, limiting influence of the separate remote wall rarefaction. These controls measure contamination at this resolution and interior window/time, not for arbitrary longer simulations. They use adaptive time steps, so their differences also include altered global step sequences.", "",
              "CFL 0.1 is compared with 0.2; CFL 0.05 is compared with 0.1. Boundary controls and time refinements must be assessed against the much larger spatial profile errors. No positivity floor or energy repair is used.", "",
              "## Work, numerical heat and solver comparison", "",
              "| C=1 radial cells | Numerical magnetic heat [J] | Relative response-energy residual | Max midpoint iterations | Advance CPU wall time [s] |", "|---:|---:|---:|---:|---:|"]
    for n in (24, 48, 96, 192):
        row = rows[f"diffusive_{n}"]; d = row["diagnostics"]
        lines.append(f"| {n} | {d['numerical_magnetic_heat_j']:.7g} | {d['response_energy_balance_relative']:.4g} | {d['max_midpoint_iterations']} | {row['advance_seconds']:.4g} |")
    coupled = [r["diagnostics"] for r in rows.values() if r["options"].get("coupled", True)]
    speed = metrics["comparisons"]["direct_vs_picard"]["advance_runtime_ratio_direct_over_picard"]
    lines += ["", f"Across coupled cases, maximum relative response-energy residual is {max(abs(d['response_energy_balance_relative']) for d in coupled):.4g}; maximum relative mass change is {max(abs(d['mass_relative_change']) for d in coupled):.4g}; maximum face divergence is {max(d['max_face_divergence_t_m'] for d in coupled):.4g} T/m. Response energy includes the full-space induced self-field, not external-field self/mutual energy or coil-circuit energy.", "",
              f"At nr=24, the retained direct method takes **{speed:.3g} times** the Picard advance time for the same C=1 fixture. The comparison table measures their field/gas profile agreement at tolerance 1e-13. Separate tests compare both exterior choices with physical and numerical diffusion independently enabled, and check source work and heat partition. These are CPU measurements from the recorded environment; setup time is stored separately. Earlier reports retain their historical source hashes and solver settings.", "",
              "![Shock profiles, convergence and numerical heating](../outputs/shock_verification/shocks.png)", "",
              "Complete per-case inputs, source/environment fingerprints, initial/final diagnostics, profiles, timings and any failure are in `outputs/shock_verification/metrics.json` and its `cases/` directory. Individual cases are saved when they finish; `--reuse` only accepts identical measurement/package provenance.", "",
              "**Remaining scope:** establish stronger shock/local-momentum accuracy, resolve coil-driven skin layers, implement a defensible evolving plasma/vacuum interface and actual end conditions, add regime-appropriate transport/ionization/Hall physics, and perform converged formation experiments. This initialized perpendicular shock is neither torus formation nor evidence of confinement or three-dimensional stability."]
    failures = [f"{r['name']} at {r['diagnostics']['time_s']:.6g} s: {r.get('failure')}" for r in rows.values() if not r["completed"]]
    if failures:
        lines += ["", "## Retained failures", ""]+["- "+failure for failure in failures]
    Path("docs/SHOCK_VERIFICATION.md").write_text("\n".join(lines)+"\n")


def verify(directory="outputs/shock_verification", *, reuse=False):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    metrics = all_measurements(directory/"cases", reuse=reuse)
    metrics["comparisons"] = comparisons(metrics)
    metrics["report_source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (directory/"metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=False)+"\n")
    plot(metrics, directory/"shocks.png"); report(metrics)
    return metrics


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse", action="store_true", help="Reuse only cases with identical package/measurement provenance")
    verify(reuse=parser.parse_args().reuse)
