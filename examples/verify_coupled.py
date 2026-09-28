"""Generate measured evidence for the fixed-chamber magnetic-vacuum coupling."""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import numpy as np
from .coupled_measurements import all_measurements


def provenance():
    digest = hashlib.sha256()
    paths = sorted(Path("plasma_bench").glob("*.py"))+[Path("examples/coupled_measurements.py"), Path(__file__)]
    for path in paths:
        digest.update(path.name.encode()); digest.update(path.read_bytes())
    return {"source_sha256": digest.hexdigest(), "python": platform.python_version(),
            "platform": platform.platform(), "versions": {k: importlib.metadata.version(k) for k in ("numpy", "scipy", "matplotlib")},
            "source_files": [str(p.relative_to(Path.cwd())) if p.is_absolute() else str(p) for p in paths],
            "initial_gas": {"density_kg_m3": 1., "pressure_pa": 1., "velocity_m_s": [0., 0.]},
            "initial_gas_exception": "passive advection has vz=0.5 m/s",
            "method": "Strang split gas SSP RK2 / adjoint magnetic implicit midpoint; exact elliptic induced vacuum",
            "boundaries": "fixed gas-filled cylinder, reflecting slip walls/ends; magnetic exterior; zero induced outer flux",
            "interpretation": "verification fixtures only; no free-interface, physical transport, formation or stability validation"}


def plot(metrics, arrays, directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    from plasma_bench.visualization import STYLE
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
        for key, error, label, color in (("force", "relative_volume_L2", "J × B force", "#ffcc75"),
                                        ("advection", "relative_flux_L2", "Flux advection", "#72ddc4"),
                                        ("diffusion", "relative_flux_L2", "Resistive diffusion", "#aa9cfc")):
            rows = metrics[key]
            axes[0, 0].loglog([r["radial_cells"] for r in rows], [r[error] for r in rows], "o-", label=label, color=color)
        axes[0, 0].set(title="Independent continuum solutions", xlabel="Radial cells", ylabel="Relative L2 error (norms differ)")
        rows = [r for r in metrics["vacuum"]["cases"] if not r["exact_boundary"]]
        axes[0, 1].loglog([r["outer_radius_and_half_length_m"] for r in rows], [r["relative_probe_L2"] for r in rows], "o-", color="#72ddc4", label="Fixed Δr=Δz=1/32 m")
        axes[0, 1].set_xticks([2, 4, 8], labels=["2", "4", "8"])
        axes[0, 1].xaxis.set_minor_formatter(NullFormatter())
        axes[0, 1].set(title="Outer boundary vs infinite-space loop integral", xlabel="Field-domain radius and half-length [m]", ylabel="Probe flux relative L2 error")
        pulse = next(r for r in metrics["pulse"] if r["radial_cells"] == 32)
        history = pulse["history"]
        time = [r["time_s"] for r in history]
        initial = history[0]["response_energy_j"]
        for key, label, color in (("source_work_j", "Applied-field work", "#72ddc4"), ("ohmic_heat_j", "Ohmic heat", "#ffcc75"),
                                  ("kinetic_energy_j", "Kinetic energy", "#aa9cfc"), ("induced_field_energy_j", "Induced self-field", "#f687a6")):
            axes[1, 0].plot(time, [r[key] for r in history], label=label, color=color)
        axes[1, 0].plot(time, [r["response_energy_j"]-initial for r in history], "--", color="white", label="Gas + induced energy gain")
        axes[1, 0].set(title="Prescribed two-coil pulse · 32 × 64 fluid cells", xlabel="Time [s]", ylabel="Integrated energy/work [J]")
        rows = [r for r in metrics["pulse"] if r["extent"] == 2 and r["configuration"]["parameters"]["cfl"] == .2]
        axes[1, 1].plot([r["radial_cells"] for r in rows], [r["diagnostics"]["maximum_speed_m_s"] for r in rows], "o-", color="#ffcc75", label="Maximum final speed")
        axes[1, 1].set(title="Pulse flow remains resolution-sensitive", xlabel="Radial cells", ylabel="Maximum speed [m/s]")
        for ax in axes.flat:
            ax.grid(alpha=.2); ax.legend(fontsize=8, labelcolor="#c9d5e5")
        fig.suptitle("DYNAMIC PLASMA + MAGNETIC VACUUM\nSmooth verification passes · fixed chamber · formation accuracy remains open", fontsize=15, weight="bold")
        fig.savefig(directory/"convergence.png", dpi=160); plt.close(fig)

        fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), layout="constrained")
        gas = arrays["gas"]; rho = gas[..., 0]
        vr, vz = gas[..., 1]/rho, gas[..., 2]/rho
        pressure = (5/3-1)*(gas[..., 3]-.5*rho*(vr*vr+vz*vz))
        for ax, value, title, units in ((axes[0], rho, "Density", "kg/m³"), (axes[1], pressure, "Thermal pressure", "Pa")):
            view = ax.pcolormesh(arrays["fluid_z_faces"], arrays["fluid_r_faces"], value, cmap="viridis", shading="flat")
            fig.colorbar(view, ax=ax, label=units, shrink=.8)
            ax.set(title=title, xlabel="z [m]", ylabel="r [m]", aspect="equal")
        stride = 3
        r = .5*(arrays["fluid_r_faces"][1:]+arrays["fluid_r_faces"][:-1])
        z = .5*(arrays["fluid_z_faces"][1:]+arrays["fluid_z_faces"][:-1])
        axes[0].quiver(z[::stride], r[::stride], vz[::stride, ::stride], vr[::stride, ::stride], color="white", alpha=.65)
        vmax = float(np.nanmax(abs(arrays["current"])))
        view = axes[2].pcolormesh(arrays["field_z"][1:-1], arrays["field_r"][1:-1], arrays["current"],
                                 cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="nearest")
        axes[2].contour(arrays["field_z"], arrays["field_r"], arrays["total_psi"], levels=18, colors="#c9d5e5", linewidths=.5)
        axes[2].plot([-1, -1, 1, 1], [0, 1, 1, 0], "--", color="white", linewidth=1)
        axes[2].set(title="Induced Jθ + total flux contours", xlabel="z [m]", ylabel="r [m]", aspect="equal", xlim=(-1.5, 1.5), ylim=(0, 1.5))
        fig.colorbar(view, ax=axes[2], label="A/m²", shrink=.8)
        fig.suptitle("FIXED-CHAMBER COIL RESPONSE · t = 0.1 s\nExterior has no fluid arrays; uncolored current region is magnetic vacuum", fontsize=14, weight="bold")
        fig.savefig(directory/"pulse.png", dpi=160); plt.close(fig)


def report(metrics):
    lines = ["# Dynamic plasma–vacuum coupling verification", "", "Generated by `python -m examples.verify_coupled`.", "",
             "This is a **fixed, gas-filled chamber reference**, not a demonstrated FRC formation model. Gas has reflecting walls/ends; the exterior contains only magnetic unknowns. External coils are prescribed analytic filaments. See [the derivation and limitations](COUPLED_NUMERICS.md).", "",
             "## Independent spatial and temporal checks", "",
             "Force uses analytic Jtheta cross B of a compact polynomial flux, integrated over actual cylindrical cells. The norm includes axis cells. Weak-field advection compares a translated compact flux at t=0.05 s in r<0.8, |z|<0.8 m, before end disturbances reach it. Diffusion compares the unbounded r² Gaussian solution at t=0.03 s in r<0.7, |z|<0.7 m; its small remote tail differs from the harmonic chamber exterior. Both passive fields have amplitude 1e-6, making fluid feedback asymptotically negligible.", "",
             "| Check | Radial cells | Relative L2 error |", "|---|---:|---:|"]
    for key, error in (("force", "relative_volume_L2"), ("advection", "relative_flux_L2"), ("diffusion", "relative_flux_L2")):
        for row in metrics[key]:
            lines.append(f"| {key} | {row['radial_cells']} | {row[error]:.7g} |")
    temporal = metrics["temporal"]
    lines += ["", f"At fixed 8×16 fluid resolution, Strang/midpoint time refinement at t=0.02 s gives momentum errors {', '.join(f'{v:.7g}' for v in temporal['momentum_errors'])} for {temporal['step_counts']} steps relative to a 256-step reference. Halving dt reduces error by factors {', '.join(f'{v:.5g}' for v in temporal['ratios'])}, consistent with second order in this smooth case.", "",
              "## Separate magnetic-domain error", "",
              "An independently integrated circular-loop Green function gives infinite-space flux from a smooth compact toroidal current. Three nonsingular probe locations are fixed. Quadrature orders 32 and 48 agree to relative " + f"{metrics['vacuum']['quadrature_32_vs_48_relative']:.3g}. Supplying that integral on the finite boundary isolates spatial error; the actual dynamic solver instead uses zero induced outer flux.", "",
              "| Cells/m | Outer radius and half-length [m] | Boundary | Probe relative L2 error |", "|---:|---:|---|---:|"]
    for row in metrics["vacuum"]["cases"]:
        lines.append(f"| {row['cells_per_m']} | {row['outer_radius_and_half_length_m']} | {'Exact infinite-space data' if row['exact_boundary'] else 'Zero induced flux'} | {row['relative_probe_L2']:.7g} |")
    lines += ["", "Zero induced current in exterior cells establishes the discrete vacuum constraint, not an infinite-space boundary. These domain errors must be checked for every proposed experiment.", "",
              "## Full prescribed-coil response", "",
              "The unit cylinder begins at rho=1 kg/m³, p=1 Pa, rest and uniform B=sqrt(mu0) T. Two radius-1.23 m loops at z=±0.37 m apply a unipolar 0.04 s rise, 0.02 s hold, 0.04 s fall. Vacuum peak center field is 0.3 times bias. eta=0.1 mu0 ohm m, gamma=5/3. A 0.002 s source time cap supplements the CFL/diffusion bound. These are numerical SI fixtures, not physical argon transport or hardware requirements.", "",
              "| Radial cells | Outer extent [m] | CFL | Applied-field work [J] | Ohmic heat [J] | Final kinetic energy [J] | Maximum speed [m/s] |", "|---:|---:|---:|---:|---:|---:|---:|"]
    for row in metrics["pulse"]:
        d = row["diagnostics"]
        lines.append(f"| {row['radial_cells']} | {row['extent']} | {row['configuration']['parameters']['cfl']:g} | {d['source_work_j']:.7g} | {d['ohmic_heat_j']:.7g} | {d['kinetic_energy_j']:.7g} | {d['maximum_speed_m_s']:.7g} |")
    diagnostics = [r['diagnostics'] for r in metrics['pulse']+metrics['advection']+metrics['diffusion']]+temporal['diagnostics']
    lines += ["", f"Maximum relative response-energy budget residual across these dynamic cases: {max(abs(d['response_energy_balance_relative']) for d in diagnostics):.3g}. Maximum relative mass change: {max(abs(d['mass_relative_change']) for d in diagnostics):.3g}. Maximum midpoint relative residual: {max(d['max_midpoint_relative_residual'] for d in diagnostics):.3g}.", "",
              "The energy budget is gas energy plus induced self-field energy, balanced against the prescribed field's work on plasma current. It excludes external-field self/mutual energy and is not a coil-circuit delivered-energy estimate. Ohmic heating is nonnegative and separately measured; kinetic work is paired with magnetic induction. No exterior fluid, cell floors or energy repair is used.", "",
              "**The pulse flow remains mesh-sensitive, and the outer boundary still materially affects the pulse.** At 16 radial cells, increasing outer extent from 2 to 4 m reduces final applied-field work by about 19%; extent 3 to 4 m still changes it by about 3.8%. This study verifies coupling and budgets, but does not establish converged compression, magnetic shocks, an evolving free plasma surface, open ends, realistic transport, Hall validity or FRC formation. The bounded reference solver's separate cold-flow and equilibrium accuracy limitations also remain open. No stability or confinement claim follows.", "",
              "![Convergence and work budget](../outputs/coupled_verification/convergence.png)", "", "![Fixed-chamber response](../outputs/coupled_verification/pulse.png)", "",
              "Complete configurations, histories, diagnostics, dependency versions and source fingerprint are saved in `outputs/coupled_verification/metrics.json`; distinct fluid/field arrays are in the locally generated `fields.npz`."]
    Path("docs/COUPLED_VERIFICATION.md").write_text("\n".join(lines)+"\n")


def verify(directory="outputs/coupled_verification"):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    metrics, arrays = all_measurements()
    metrics["provenance"] = provenance()
    (directory/"metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=False)+"\n")
    np.savez_compressed(directory/"fields.npz", **arrays)
    plot(metrics, arrays, directory)
    report(metrics)
    return metrics


if __name__ == "__main__":
    verify()
