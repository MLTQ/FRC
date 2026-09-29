"""Scientific figures for boundary accuracy, keeping different evidence separate."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import numpy as np
from .boundary_analysis import read_case


def plot(metrics, directory):
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.6), layout="constrained")
    colors = {"standard": "#2767ad", "strip": "#c44374", "physical": "#4c8655"}
    labels = {"standard": "Second-order boundary", "strip": "First-order boundary strip", "physical": "Physical diffusion only"}
    ax = axes[0, 0]
    for method in ("standard", "strip"):
        rows = metrics["shocks"][method]
        ax.loglog([r["radial_cells"] for r in rows],
                  [abs(r["diagnostics"]["response_energy_balance_residual_j"]) for r in rows],
                  "o-", color=colors[method], label=labels[method])
    ax.set(title="A  Original shock: response energy", xlabel="Radial cells", ylabel="Absolute residual [J]")
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    for method in ("physical", "standard", "strip"):
        rows = [r for r in metrics["skin"] if r["method"] == method and r["cfl"] == .2]
        ax.loglog([r["radial_cells"] for r in rows], [r["maximum_B_error"] for r in rows],
                  "o-", color=colors[method], label=labels[method])
    ax.set(title="B  Independent cylindrical skin solution", xlabel="Radial cells", ylabel="Maximum B L2 / applied amplitude")
    ax.legend(fontsize=8)
    ax = axes[0, 2]
    for method in ("standard", "strip"):
        rows = [r for r in metrics["pulse_summary"] if r["name"] in [f"{method}_{n}" for n in (16, 32, 64, 128)]]
        ax.plot([r["radial_cells"] for r in rows], [r["samples"][-1]["maximum_speed_m_s"] for r in rows],
                "o-", color=colors[method], label=labels[method])
    ax.set(title="C  Same coil pulse: final peak speed", xlabel="Radial cells", ylabel="Peak speed [m/s]")
    ax.legend(fontsize=8)
    ax = axes[1, 0]
    for method in ("standard", "strip"):
        rows = [r for r in metrics["pulse_refinement"] if r["coarse_or_test"].startswith(method+"_")]
        ax.semilogx([int(r["coarse_or_test"].split("_")[-1]) for r in rows],
                  [100*max(t["velocity_relative_L2"] for t in r["history"]) for r in rows],
                  "o-", color=colors[method], label=labels[method])
    ax.set(title="D  Pulse flow: adjacent-grid difference", xlabel="Coarse radial cells (fine = 2×)", ylabel="Maximum velocity L2 difference [%]")
    ax.legend(fontsize=8)
    ax = axes[1, 1]
    for method in ("standard", "strip"):
        for n, style, alpha in ((32, ":", .7), (64, "--", .85), (128, "-", 1.)):
            row, fields = read_case(directory/"pulses", f"{method}_{n}")
            current = fields["t8_current"]
            nz = row["fluid_mesh"]["nz"]
            z = np.linspace(-1, 1, nz+1)
            profile = np.mean(current[:, abs(z) <= .1], axis=1)
            ax.plot(np.linspace(0, 1, n+1), profile, style, color=colors[method], alpha=alpha,
                    label=f"{method}, {n}")
    ax.set(title="E  Pulse current at end of rise", xlabel="Radius [m]", ylabel="Central-band Jθ [A/m²]", xlim=(.65, 1))
    ax.legend(fontsize=7, ncol=2)
    ax = axes[1, 2]
    for method in ("standard", "strip"):
        rows = [r for r in metrics["pulse_summary"] if r["name"] in [f"{method}_{n}" for n in (16, 32, 64, 128)]]
        ax.loglog([r["radial_cells"] for r in rows],
                  [100*r["diagnostics"]["numerical_interface_heat_j"]/r["diagnostics"]["ohmic_heat_j"] for r in rows],
                  "o-", color=colors[method], label=labels[method])
    ax.set(title="F  Pulse boundary dissipation", xlabel="Radial cells", ylabel="Numerical interface / physical heat [%]")
    ax.legend(fontsize=8)
    for ax, ticks in ((axes[0, 0], [24, 48, 96, 192, 384]),
                      (axes[0, 1], [16, 32, 64, 128]),
                      (axes[0, 2], [16, 32, 64, 128]),
                      (axes[1, 0], [16, 32, 64]),
                      (axes[1, 2], [16, 32, 64, 128])):
        ax.set_xticks(ticks, [str(n) for n in ticks])
        ax.xaxis.set_minor_locator(NullLocator())
    axes[1, 2].set_yticks([.5, 1, 2, 5, 10, 20], ["0.5", "1", "2", "5", "10", "20"])
    axes[1, 2].yaxis.set_minor_locator(NullLocator())
    for ax in axes.flat:
        ax.grid(alpha=.2)
    fig.suptitle("Boundary verification: shock, analytic skin response and driven gas\nFixed-chamber / kinematic references; no formation or free-surface claim", fontsize=14)
    fig.savefig(directory/"boundary_verification.png", dpi=170)
    plt.close(fig)
