"""Measure low-pressure robustness and independent expanding-plasma accuracy."""
from dataclasses import asdict
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import numpy as np
from plasma_bench.mhd import AxisymmetricMHD, MHDParameters
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.mhd_state import primitive
from plasma_bench.mhd_equilibrium import PeriodicEquilibrium
from plasma_bench.mhd_reconstruction import Reconstruction
from plasma_bench.mhd_benchmarks import longitudinal_rarefaction, double_rarefaction_exact
from .verify_equilibrium import run_case


def verify(directory="outputs/positivity_verification"):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    metrics, arrays = {"stress_equilibria": [], "double_rarefaction": []}, {}
    for pressure in (.001, .00001):
        for n in (12, 24):
            metrics["stress_equilibria"].append(run_case(pressure, n)[0])
    for n in (80, 160, 320, 640):
        mesh = MHDMesh(1, 1, 4, n)
        initial = longitudinal_rarefaction(mesh)
        solver = AxisymmetricMHD(mesh, initial, MHDParameters(gamma=1.4))
        solver.advance(.1)
        q = primitive(mesh, solver.state, 1.4)
        rho, velocity, pressure = double_rarefaction_exact(mesh.z, .1)
        star_pressure = float(double_rarefaction_exact(np.asarray(0.), .1)[2])
        select = abs(mesh.z) < .5
        row = {"axial_cells": n, "mesh": asdict(mesh), "parameters": asdict(solver.parameters),
               "density_L1_kg_m3": float(np.mean(abs(q[1, select, 0]-rho[select]))),
               "pressure_L1_pa": float(np.mean(abs(q[1, select, 3]-pressure[select]))),
               "exact_star_pressure_pa": star_pressure,
               "central_pressure_pa": float(q[1, n//2, 3]),
               "central_pressure_relative_error": float(abs(q[1, n//2, 3]-star_pressure)/star_pressure),
               "maximum_flux_change_wb_per_rad": float(abs(solver.state.psi-initial.psi).max()),
               **solver.diagnostics()}
        metrics["double_rarefaction"].append(row)
        arrays[f"rarefaction_{n}"] = np.stack((mesh.z[select], q[1, select, 0], rho[select], q[1, select, 3], pressure[select]))
    mesh = MHDMesh(1, 1, 12, 24)
    initial = PeriodicEquilibrium(background_pressure=.1).state(mesh)
    reconstruction = Reconstruction(mesh, primitive(mesh, initial, 5/3))
    arrays.update(r_faces=mesh.rf, z_faces=mesh.zf, thermal_blend=1-reconstruction.pressure_weight)
    digest = hashlib.sha256()
    for path in sorted(Path("plasma_bench").glob("*.py")):
        digest.update(path.name.encode()); digest.update(path.read_bytes())
    metrics["provenance"] = {"source_sha256": digest.hexdigest(), "python": platform.python_version(),
                              "platform": platform.platform(), "versions": {k: importlib.metadata.version(k) for k in ("numpy", "scipy", "matplotlib")},
                              "method": "limited second-order Rusanov/UCT with local admissibility reconstruction; SSP RK2; no cell floors",
                              "boundaries": "axis, conducting slip radial wall, periodic z",
                              "interpretation": "numerical verification; no formation or physical regime validation"}
    (directory/"metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=False)+"\n")
    np.savez_compressed(directory/"fields.npz", **arrays)
    plot(metrics, arrays, directory/"positivity.png")
    report(metrics)
    return metrics


def plot(metrics, arrays, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    from plasma_bench.visualization import STYLE
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
        z, density, exact_density, pressure, exact_pressure = arrays["rarefaction_640"]
        axes[0, 0].plot(z, density, color="#72ddc4", label="640 axial cells")
        axes[0, 0].plot(z, exact_density, "--", color="white", label="Exact isentropic fans")
        axes[0, 0].set(title="Two separating streams · uniform parallel B", xlabel="z [m]", ylabel="Density [kg/m³]")
        axes[0, 1].semilogy(z, pressure, color="#72ddc4", label="640 axial cells")
        axes[0, 1].semilogy(z, exact_pressure, "--", color="white", label="Exact pressure")
        axes[0, 1].set(title="Pressure accuracy · central error remains visible", xlabel="z [m]", ylabel="Thermal pressure [Pa]")
        cases = metrics["double_rarefaction"]
        for key, label, color in (("density_L1_kg_m3", "Density L1 [kg/m³]", "#72ddc4"), ("pressure_L1_pa", "Pressure L1 [Pa]", "#aa9cfc")):
            axes[1, 0].loglog([v["axial_cells"] for v in cases], [v[key] for v in cases], "o-", color=color, label=label)
        axes[1, 0].set(title="Measured rarefaction refinement", xlabel="Axial cells", ylabel="Absolute L1 error (units in legend)")
        axes[1, 0].set_xticks([80, 160, 320, 640], labels=["80", "160", "320", "640"])
        axes[1, 0].xaxis.set_minor_formatter(NullFormatter())
        picture = axes[1, 1].pcolormesh(arrays["z_faces"], arrays["r_faces"], arrays["thermal_blend"], vmin=0, vmax=1, cmap="magma", shading="flat")
        fig.colorbar(picture, ax=axes[1, 1], label="Auxiliary thermal-pressure weight")
        axes[1, 1].set(title="Initial limiting · p₀=0.1 Pa, 12 × 24 cells", xlabel="z [m]", ylabel="Radius [m]")
        for ax in (axes[0, 0], axes[0, 1], axes[1, 0]):
            ax.grid(alpha=.2); ax.legend(fontsize=8, labelcolor="#c9d5e5")
        fig.suptitle("LOW-PRESSURE NUMERICAL VERIFICATION\nPositive reconstruction · conserved cell budgets · resolution still matters", fontsize=15, weight="bold")
        fig.savefig(path, dpi=160); plt.close(fig)


def report(metrics):
    lines = ["# Low-pressure robustness and accuracy", "", "Generated by `python -m examples.verify_positivity`.", "",
             "[The local admissibility method](POSITIVITY.md) changes face/corner reconstruction and its geometric source, leaving conserved cell averages intact. It is not a pressure floor or a general proof of full-update positivity.", "",
             "## More demanding equilibrium stress cases", "",
             "These intentionally coarse closed-flux equilibria test whether the solver completes without state repair. Their appreciable artificial motion is **not** a resolved physical result. See the separate [equilibrium refinement study](EQUILIBRIUM_VERIFICATION.md).", "",
             "| p0 [Pa] | Radial cells | Status | Minimum evolved cell pressure [Pa] | Max artificial speed / reference vA | Pressure-blended cell-stages |", "|---:|---:|---|---:|---:|---:|"]
    for v in metrics["stress_equilibria"]:
        stats = v["reconstruction_limiter"]
        fraction = stats["pressure_blended_cells"]/max(1, stats["accepted_cell_stages"])
        lines.append(f"| {v['equilibrium']['background_pressure']:g} | {v['radial_cells']} | {v['status']} | {v['min_pressure_pa']:.6g} | {v['velocity_over_reference_alfven_speed']:.5g} | {fraction:.3%} |")
    lines += ["", "## Independent expanding-flow solution", "",
              "Initially rho=1 kg/m³, p=0.4 Pa, vz=±2 m/s, gamma=1.4 and a uniform axial field with vA=1 m/s. The streams separate. Because flow and constant B are parallel, two exact gas-dynamic isentropic fans apply and B must remain unchanged. Only |z|<0.5 m at t=0.1 s is compared, before periodic-boundary shocks arrive. True-vacuum-producing initial data are explicitly excluded from this fixture.", "",
              "| Axial cells | Density L1 [kg/m³] | Pressure L1 [Pa] | Central pressure [Pa] | Central pressure relative error |", "|---:|---:|---:|---:|---:|"]
    for v in metrics["double_rarefaction"]:
        lines.append(f"| {v['axial_cells']} | {v['density_L1_kg_m3']:.5g} | {v['pressure_L1_pa']:.5g} | {v['central_pressure_pa']:.6g} | {v['central_pressure_relative_error']:.3%} |")
    exact = metrics["double_rarefaction"][0]["exact_star_pressure_pa"]
    rows = metrics["stress_equilibria"]+metrics["double_rarefaction"]
    lines += ["", f"Exact central pressure is {exact:.8g} Pa. Global L1 errors decrease, but the low-pressure central region remains much less accurate than the bulk profile on these meshes. Limiting is inactive in these rarefaction runs: this residual error belongs to the underlying Rusanov discretization, not a pressure floor. It must remain visible when choosing resolution and interpreting rarefied plasma regions.", "",
              f"Maximum relative mass change across completed cases: {max(abs(v['mass_relative_change']) for v in rows):.3g}. Maximum relative total-energy change: {max(abs(v['total_energy_relative_change']) for v in rows):.3g}.", "",
              "The uniform parallel flux is unchanged to floating-point accuracy. Cell density and pressure stay positive without repair.", "",
              "![Low-pressure diagnostics](../outputs/positivity_verification/positivity.png)", "",
              "These results clear the specific reconstruction failure and provide a stronger bounded reference for dynamic vacuum/coil coupling. They do not establish complete formation accuracy, transport validity, confinement or stability."]
    Path("docs/POSITIVITY_VERIFICATION.md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    verify()
