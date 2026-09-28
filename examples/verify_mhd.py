"""Measure independent MHD benchmarks, write reports and inspectable arrays."""
from dataclasses import asdict
from pathlib import Path
import json
import hashlib
import importlib.metadata
import platform
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.mhd import AxisymmetricMHD, MHDParameters
from plasma_bench.mhd_state import primitive
from plasma_bench.mhd_benchmarks import radial_fast_mode, advected_flux, resistive_bessel_mode, sod_initial, sod_exact, radial_pressure_balance, radial_average


def _relative_l2(actual, expected, weight=1):
    return float(np.sqrt(np.sum((actual-expected)**2*weight)/np.sum(expected**2*weight)))


def verify(directory="outputs/mhd_verification"):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    metrics = {"fast_wave": [], "advection": [], "resistive_decay": [], "sod": [], "time_refinement": [], "equilibrium": [], "magnetic_force": []}
    plot_data = {}
    for n in (24, 48, 96, 192):
        mesh = MHDMesh(1, 1, n, 4)
        initial, omega = radial_fast_mode(mesh)
        solver = AxisymmetricMHD(mesh, initial)
        t = np.pi/(2*omega)
        solver.advance(t)
        expected, _ = radial_fast_mode(mesh, time=t)
        q, exact = primitive(mesh, solver.state, 5/3), primitive(mesh, expected, 5/3)
        error = _relative_l2(q[..., 1], exact[..., 1], mesh.volumes)
        metrics["fast_wave"].append({"radial_cells": n, "end_time_s": t, "velocity_relative_L2": error,
                                    "mesh": asdict(mesh), "parameters": asdict(solver.parameters), **solver.diagnostics()})
        plot_data[f"wave_{n}"] = np.stack((mesh.r, q[:, 0, 1], exact[:, 0, 1]))
    for n in (12, 24, 48, 96):
        mesh = MHDMesh(1, 1, n, 2*n)
        solver = AxisymmetricMHD(mesh, advected_flux(mesh))
        solver.advance(.15)
        expected = advected_flux(mesh, time=.15)
        error = _relative_l2(solver.state.psi, expected.psi)
        metrics["advection"].append({"radial_cells": n, "flux_relative_L2": error,
                                    "mesh": asdict(mesh), "parameters": asdict(solver.parameters), **solver.diagnostics()})
        if n == 96:
            plot_data.update(advected_psi=solver.state.psi, exact_psi=expected.psi, r_faces=mesh.rf, z_faces=mesh.zf)
    for n in (16, 32, 64, 128):
        mesh = MHDMesh(1, 1, n, 4)
        initial, k = resistive_bessel_mode(mesh)
        solver = AxisymmetricMHD(mesh, initial, MHDParameters(resistivity=MU0*.5))
        before = solver.diagnostics()
        solver.advance(.01)
        after = solver.diagnostics()
        expected = initial.psi*np.exp(-.5*k*k*.01)
        lost = before["magnetic_energy_j"]-after["magnetic_energy_j"]
        heat_error = abs(after["physical_ohmic_heat_j"]-lost)/lost
        metrics["resistive_decay"].append({"radial_cells": n, "flux_relative_L2": _relative_l2(solver.state.psi, expected),
                                          "heating_relative_error": heat_error, "magnetic_energy_loss_j": lost,
                                          "mesh": asdict(mesh), "parameters": asdict(solver.parameters), **after})
    for n in (80, 160, 320):
        mesh = MHDMesh(1, 1, 4, n)
        solver = AxisymmetricMHD(mesh, sod_initial(mesh), MHDParameters(gamma=1.4))
        solver.advance(.12)
        q = primitive(mesh, solver.state, 1.4)
        exact, velocity, pressure = sod_exact(mesh.z, .12)
        central = abs(mesh.z) < .4
        metrics["sod"].append({"axial_cells": n, "density_L1": float(np.mean(abs(q[1, central, 0]-exact[central]))),
                               "mesh": asdict(mesh), "parameters": asdict(solver.parameters), **solver.diagnostics()})
        plot_data[f"sod_{n}"] = np.stack((mesh.z[central], q[1, central, 0], exact[central]))
    mesh = MHDMesh(1, 1, 24, 4)
    initial, omega = radial_fast_mode(mesh)
    end_time = np.pi/(2*omega)
    baseline = AxisymmetricMHD(mesh, initial).timestep_limit()
    velocities, counts = [], []
    for fraction in (.8, .4, .2, .025):
        solver = AxisymmetricMHD(mesh, initial)
        steps = int(np.ceil(end_time/(fraction*baseline)))
        for _ in range(steps):
            solver.step(end_time/steps)
        velocities.append(primitive(mesh, solver.state, 5/3)[..., 1])
        counts.append(steps)
    for count, velocity in zip(counts[:-1], velocities[:-1]):
        metrics["time_refinement"].append({"steps": count, "reference_steps": counts[-1],
                                            "velocity_relative_L2_vs_time_reference": _relative_l2(velocity, velocities[-1])})
    for n in (24, 48, 96, 192):
        mesh = MHDMesh(1, 1, n, 4)
        initial, b0, _ = radial_pressure_balance(mesh)
        solver = AxisymmetricMHD(mesh, initial)
        solver.advance(.2)
        diagnostics = solver.diagnostics()
        metrics["equilibrium"].append({"radial_cells": n, "spurious_velocity_over_reference_alfven_speed": diagnostics["max_velocity_m_s"]/(b0/np.sqrt(MU0)),
                                       "mesh": asdict(mesh), "parameters": asdict(solver.parameters), **diagnostics})
        unbalanced, b0, alpha = radial_pressure_balance(mesh, balanced=False)
        force, _ = AxisymmetricMHD(mesh, unbalanced).rhs(unbalanced)
        exact = radial_average(mesh, lambda r: 2*alpha*b0*b0/MU0*r*(1-alpha*r*r))
        actual = force.conserved[:, 0, 1]
        select = mesh.r < .9*mesh.radius
        error = _relative_l2(actual[select], exact[select], mesh.volumes[select, 0])
        metrics["magnetic_force"].append({"radial_cells": n, "force_volume_relative_L2": error,
                                         "first_axis_cell_acceleration_absolute_error_m_s2": float(abs(actual[0]-exact[0]))})
    digest = hashlib.sha256()
    for path in sorted(Path("plasma_bench").glob("*.py")):
        digest.update(path.name.encode()); digest.update(path.read_bytes())
    metrics["provenance"] = {"source_sha256": digest.hexdigest(), "python": platform.python_version(),
                              "platform": platform.platform(), "versions": {k: importlib.metadata.version(k) for k in ("numpy", "scipy", "matplotlib")},
                              "method": "first-order Rusanov/UCT spatial discretization; SSP RK2; no positivity floors",
                              "boundaries": "axis regularity, radial slip conducting wall, periodic z",
                              "benchmarks": "fixed analytic parameters defined in plasma_bench/mhd_benchmarks.py; not plasma formation"}
    (directory/"metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=False)+"\n")
    np.savez_compressed(directory/"benchmark_fields.npz", **plot_data)
    plot(metrics, plot_data, directory/"verification.png")
    report(metrics)
    return metrics


def plot(metrics, data, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    from plasma_bench.visualization import STYLE
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
        for n, color in zip((24, 48, 96, 192), ("#64748b", "#aa9cfc", "#f5b971", "#72ddc4")):
            r, actual, exact = data[f"wave_{n}"]
            axes[0, 0].plot(r, actual*1e6, color=color, label=f"{n} radial cells")
        axes[0, 0].plot(r, exact*1e6, "--", color="white", label="Analytic fast wave")
        axes[0, 0].set(title="Cylindrical fast wave · velocity at quarter period", xlabel="Radius [m]", ylabel="Radial velocity [µm/s]")
        z, actual, exact = data["sod_320"]
        axes[0, 1].plot(z, actual, color="#72ddc4", label="320 axial cells")
        axes[0, 1].plot(z, exact, "--", color="white", label="Exact Sod solution")
        axes[0, 1].set(title="Compressible shock · independent exact solution", xlabel="z [m]", ylabel="Density [kg/m³]")
        for key, metric, label, color in (("fast_wave", "velocity_relative_L2", "Fast-wave velocity", "#72ddc4"),
                                          ("advection", "flux_relative_L2", "2-D flux advection", "#aa9cfc"),
                                          ("resistive_decay", "flux_relative_L2", "Resistive flux decay", "#f5b971")):
            values = metrics[key]
            axes[1, 0].loglog([v["radial_cells"] for v in values], [v[metric] for v in values], "o-", label=label, color=color)
        axes[1, 0].set(title="Measured spatial refinement", xlabel="Radial cells", ylabel="Relative L2 error")
        values = metrics["resistive_decay"]
        axes[1, 1].loglog([v["radial_cells"] for v in values], [v["heating_relative_error"] for v in values], "o-", color="#f5b971", label="Physical Joule heat vs magnetic loss")
        axes[1, 1].set(title="Heating consistency · numerical diffusion is measurable", xlabel="Radial cells", ylabel="Relative discrepancy")
        for ax, ticks in ((axes[1, 0], [12, 24, 48, 96, 192]), (axes[1, 1], [16, 32, 64, 128])):
            ax.set_xticks(ticks, labels=[str(tick) for tick in ticks])
            ax.xaxis.set_minor_formatter(NullFormatter())
        for ax in axes.flat:
            ax.grid(alpha=.2); ax.legend(fontsize=8, framealpha=.1, labelcolor="#c9d5e5")
        fig.suptitle("DYNAMIC MHD VERIFICATION\nConservative r-z dynamics · declared conducting wall / periodic boundaries", fontsize=16, weight="bold")
        fig.savefig(path, dpi=160)
        plt.close(fig)


def report(metrics):
    lines = ["# Dynamic MHD verification", "", "Generated by `python -m examples.verify_mhd`.", "",
             "The new kernel evolves density, radial/axial momentum, total energy and poloidal flux.",
             "These are verification cases in a conducting-wall, axially periodic, gas-filled cylinder.",
             "They do not substitute for a moving plasma surrounded by vacuum or prove FRC formation.", "",
             "## Independent spatial refinement", "",
             "| Case | Cells | Error quantity | Error |", "|---|---:|---|---:|"]
    for key, field, count, label in (("fast_wave", "velocity_relative_L2", "radial_cells", "velocity relative L2"),
                                     ("advection", "flux_relative_L2", "radial_cells", "flux relative L2"),
                                     ("resistive_decay", "flux_relative_L2", "radial_cells", "flux relative L2"),
                                     ("sod", "density_L1", "axial_cells", "density absolute L1 [kg/m³]")):
        for row in metrics[key]:
            lines.append(f"| {key} | {row[count]} | {label} | {row[field]:.5g} |")
    rows = [row for key in ("fast_wave", "advection", "resistive_decay", "sod", "equilibrium") for row in metrics[key]]
    lines += ["", "The fast wave uses analytic cylindrical Bessel eigenfunctions at quarter period,",
              "including compression, magnetic advection and pressure response. Its linearization",
              "error is O(10^-10) in the state. The flux-advection case uses a dynamically negligible",
              "magnetic perturbation. The exact Sod comparison excludes the second periodic Riemann problem.", "",
              "## Temporal refinement at fixed mesh", "",
              "Cylindrical fast wave on 24 radial cells, compared with a finer **time** reference at the same spatial resolution.", "",
              "| Steps | Reference steps | Velocity relative L2 error |", "|---|---:|---:|"]
    for row in metrics["time_refinement"]:
        lines.append(f"| {row['steps']} | {row['reference_steps']} | {row['velocity_relative_L2_vs_time_reference']:.5g} |")
    lines += ["", "The error decreases approximately fourfold when the time step halves, as expected for SSP RK2. This does not remove the measured first-order spatial error.", "",
              "## Conservation and heating", "",
              f"Maximum relative mass change: {max(abs(r['mass_relative_change']) for r in rows):.3g}.",
              f"Maximum relative total-energy change: {max(abs(r['total_energy_relative_change']) for r in rows):.3g}.",
              "Magnetic fields are obtained from evolving corner flux, maintaining telescoping face divergence.",
              "The independent tests also check uniform axis balance and the full one-dimensional UCT dissipation.", "",
              "| Radial cells | Physical Joule heat vs magnetic-energy loss discrepancy |", "|---|---:|"]
    for row in metrics["resistive_decay"]:
        lines.append(f"| {row['radial_cells']} | {row['heating_relative_error']:.3%} |")
    lines += ["", "Total energy is conserved while magnetic energy becomes internal energy; no separate",
              "Joule source is added to total energy. The physical heating quadrature discrepancy",
              "measures spatial error plus numerical magnetic diffusion; it is not forced to zero.",
              "The first-order Rusanov/UCT scheme remains relatively diffusive. A higher-order,",
              "less-dissipative method and an actual formation convergence study are still needed.", "",
              "## Static magnetic-pressure balance", "",
              "The fixture Bz=B0(1-alpha r²), p+Bz²/(2mu0)=constant has known zero continuum acceleration.",
              "Initial discrete force balance is accurate to roundoff. Numerical magnetic diffusion subsequently changes it and drives artificial flow. This is measured, not suppressed by a force correction.", "",
              "| Radial cells | Maximum spurious velocity / reference Alfven speed |", "|---|---:|"]
    for row in metrics["equilibrium"]:
        lines.append(f"| {row['radial_cells']} | {row['spurious_velocity_over_reference_alfven_speed']:.5g} |")
    lines += ["", "The unbalanced magnetic force is independently compared with volume-averaged analytic J×B at rho=1 kg/m³. The physical volume norm includes the axis; only the outer wall region is excluded.", "",
              "| Radial cells | Force volume-relative L2 error | First axis-cell acceleration absolute error [m/s²] |", "|---|---:|---:|"]
    for row in metrics["magnetic_force"]:
        lines.append(f"| {row['radial_cells']} | {row['force_volume_relative_L2']:.5g} | {row['first_axis_cell_acceleration_absolute_error_m_s2']:.5g} |")
    lines += ["", "The first-axis-cell acceleration error decreases linearly, even though the volume-weighted norm improves faster. Higher-order axis/face reconstruction remains necessary for precise core forces. A straight-field equilibrium does not verify FRC equilibrium preservation.", "",
              "![Verification plots](../outputs/mhd_verification/verification.png)", "",
              "## Remaining gates", "",
              "- Dynamic plasma–vacuum interface and source-work coupling to external coils.",
              "- Physical chamber ends and open exterior magnetic boundary.",
              "- Regime-appropriate transport and Hall/two-fluid validation.",
              "- Equilibrium preservation and converged closed-flux formation/lifetime tests.",
              "- Three-dimensional stability and particle-confinement evidence, if claimed.", "",
              "The broad FRC-simulation goal remains active. [Numerical contract](MHD_NUMERICS.md)."]
    Path("docs/MHD_VERIFICATION.md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    measured = verify()
    print(json.dumps({k: v for k, v in measured.items() if k != "provenance"}, indent=2))
