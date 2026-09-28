"""Scientific plots and self-contained, playable time-slice viewers."""
from pathlib import Path
import html
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import ListedColormap, BoundaryNorm
from .config import Experiment


COLORS = {"unipolar": "#72ddc4", "rlc": "#f5b971", "burst": "#aa9cfc"}
STYLE = {"figure.facecolor": "#101722", "axes.facecolor": "#151f2c", "savefig.facecolor": "#101722",
         "text.color": "#eef3fa", "axes.labelcolor": "#c9d5e5", "xtick.color": "#a4b5c9",
         "ytick.color": "#a4b5c9", "axes.edgecolor": "#475569", "grid.color": "#64748b",
         "font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12}


def _slice_data(result):
    a = result.arrays
    config = Experiment.from_dict(result.metadata["configuration"])
    ri = np.flatnonzero(a["r"] <= config.chamber.radius+1e-12)
    zi = np.flatnonzero(abs(a["z"]) <= config.chamber.half_length+1e-12)
    return config, ri, zi


def _frame_figure(result):
    a = result.arrays
    config, ri, zi = _slice_data(result)
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), gridspec_kw={"height_ratios": [1.4, 1]}, layout="constrained")
    extent = [a["z"][zi[0]]*100, a["z"][zi[-1]]*100, 0, a["r"][ri[-1]]*100]
    bz = a["bz"][:, ri][:, :, zi]*1e3
    j = a["j_theta"][:, ri][:, :, zi]*1e-6
    bmax = max(float(np.abs(bz).max()), 1e-9)
    jmax = max(float(np.abs(j).max()), 1e-9)
    im0 = axes[0, 0].imshow(bz[0], origin="lower", extent=extent, aspect="auto", cmap="RdBu_r", vmin=-bmax, vmax=bmax)
    im1 = axes[0, 1].imshow(j[0], origin="lower", extent=extent, aspect="auto", cmap="PuOr_r", vmin=-jmax, vmax=jmax)
    fig.colorbar(im0, ax=axes[0, 0], label="Bz [mT]", shrink=.85)
    fig.colorbar(im1, ax=axes[0, 1], label="Jθ [MA/m²]", shrink=.85)
    for ax in axes[0]:
        ax.set(xlabel="z [cm]", ylabel="r [cm]")
        p = config.material.region
        ax.plot([-p.half_length*100, -p.half_length*100, p.half_length*100, p.half_length*100],
                [0, p.radius*100, p.radius*100, 0], "--", color="#a8fff1", lw=1, alpha=.8)
    axes[0, 0].set_title("Axial field + ψ contours · O/X markers")
    axes[0, 1].set_title("Induced toroidal current · fixed plasma outline")
    t = a["times"]*1e6
    axes[1, 0].plot(t, a["center_bz"]*1e3, color="#72ddc4", label="Center probe · total")
    axes[1, 0].plot(t, a["vacuum_center_bz"]*1e3, "--", color="#a4b5c9", label="Center probe · vacuum")
    axes[1, 0].plot(t, a["edge_bz"]*1e3, color="#f5b971", alpha=.8, label="Radial probe · total")
    axes[1, 0].axhline(0, color="white", alpha=.3, lw=.7)
    axes[1, 0].set(xlabel="Time [µs]", ylabel="Bz [mT]", title="Virtual probes")
    axes[1, 1].plot(t, a["induced_magnetic_energy"], color="#72ddc4", label="Stored induced field")
    axes[1, 1].plot(t, a["ohmic_heat"], color="#f5b971", label="Accumulated Ohmic heat")
    axes[1, 1].plot(t, a["numerical_loss"], "--", color="#aa9cfc", label="Numerical loss")
    axes[1, 1].set(xlabel="Time [µs]", ylabel="Energy [J]", title="Induced-field energy accounting")
    lines = [ax.axvline(0, color="#eef3fa", alpha=.7, lw=1) for ax in axes[1]]
    for ax in axes[1]:
        ax.grid(alpha=.18)
        ax.legend(fontsize=8, framealpha=.15, labelcolor="#c9d5e5")
    title = fig.suptitle("", fontsize=15, weight="bold")
    contour = [None]
    markers = []

    def update(i):
        im0.set_data(bz[i]); im1.set_data(j[i])
        if contour[0] is not None:
            contour[0].remove()
        for marker in markers:
            marker.remove()
        markers.clear()
        p = a["psi"][i][np.ix_(ri, zi)]
        topo = result.metadata["topology"][i]
        levels = list(np.linspace(float(p.min()), float(p.max()), 22)[1:-1])
        for point in topo.get("o_points", []):
            levels.extend(point["psi_wb"]*fraction for fraction in (.1, .5, .9))
        levels = np.unique([v for v in levels if p.min() < v < p.max()])
        if np.ptp(p) > 0:
            contour[0] = axes[0, 0].contour(a["z"][zi]*100, a["r"][ri]*100, p, levels=levels, colors="#16283a", linewidths=.65, alpha=.8)
            contour[0].set_path_effects([path_effects.withStroke(linewidth=1.2, foreground="#f8fafc")])
        for key, symbol, color in (("o_points", "o", "#8a5600"), ("x_points", "x", "#a6003d")):
            for point in topo.get(key, []):
                markers.extend(axes[0, 0].plot(point["z_m"]*100, point["r_m"]*100, marker=symbol, color=color, ms=5, fillstyle="none"))
        for line in lines:
            line.set_xdata([a["frame_times"][i]*1e6]*2)
        title.set_text(f"PLASMA BENCH  /  {config.label.upper()}  /  Stage {config.stage}  /  {a['frame_times'][i]*1e6:.1f} µs\n"
                       f"Closed flux in prescribed plasma: {'detected' if topo['closed_flux'] else 'not detected'}   ·   FRC formation not established")
        return im0, im1, title

    return fig, update


def snapshot(result, path, frame=None):
    if frame is None:
        # Show largest absolute induced current; vacuum snapshot uses largest source excursion.
        metric = np.max(np.abs(result.arrays["j_theta"]), axis=(1, 2))
        if metric.max() == 0:
            idx = result.arrays["frame_indices"]
            metric = np.abs(result.arrays["center_bz"][idx]-result.arrays["center_bz"][0])
        frame = int(metric.argmax())
    with plt.rc_context(STYLE):
        fig, update = _frame_figure(result)
        update(frame)
        fig.savefig(path, dpi=160)
        plt.close(fig)


def viewer(result, path, max_frames=51, gif=None):
    """Export a self-contained HTML movie with play, pause, step and time slider."""
    indices = np.unique(np.linspace(0, len(result.arrays["frame_times"])-1,
                                    min(max_frames, len(result.arrays["frame_times"])), dtype=int))
    with plt.rc_context({**STYLE, "animation.embed_limit": 100}):
        fig, update = _frame_figure(result)
        animation = FuncAnimation(fig, update, frames=indices, interval=100, blit=False)
        body = animation.to_jshtml(default_mode="loop")
        # Matplotlib's icon stylesheet is a CDN dependency; replace it for offline use.
        body = re.sub(r'<link[^>]*font-awesome[^>]*>', '', body, flags=re.DOTALL)
        for icon, label in (("minus", "Slower"), ("fast-backward", "First"),
                            ("step-backward", "Back"), ("pause", "Pause"),
                            ("step-forward", "Next"), ("fast-forward", "Last"), ("plus", "Faster")):
            body = body.replace(f'<i class="fa fa-{icon}"></i>', label)
        body = body.replace('<i class="fa fa-play fa-flip-horizontal"></i>', 'Reverse')
        body = body.replace('<i class="fa fa-play"></i>', 'Play')
        if gif is not None:
            animation.save(gif, writer=PillowWriter(fps=10), dpi=100)
        warnings = "".join(f"<li>{html.escape(w)}</li>" for w in result.metadata["warnings"])
        document = ("<!doctype html><html><meta charset='utf-8'><title>Plasma Bench · Time viewer</title>"
                    "<style>body{background:#101722;color:#dbe6f3;font:16px system-ui;max-width:1280px;margin:30px auto;padding:0 20px}"
                    "h1{font-size:24px}li{margin:8px 0}.animation{max-width:100%}img{max-width:100%;height:auto}button{cursor:pointer}</style>"
                    "<h1>Virtual plasma bench · verified model, exploratory parameters</h1>"
                    "<p>Play or scrub saved time slices. Fields use SI internally. The dashed boundary is the fixed plasma region.</p>"
                    +body+"<details><summary>Model assumptions and validity warnings</summary><ul>"+warnings+"</ul></details></html>")
        Path(path).write_text(document)
        plt.close(fig)


def comparison_plot(results, path):
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
        for name, result in results.items():
            a, color = result.arrays, COLORS[name]
            t = a["times"]*1e6
            axes[0, 0].plot(t, a["center_bz"]*1e3, color=color, label=name)
            axes[0, 0].plot(t, a["vacuum_center_bz"]*1e3, "--", color=color, alpha=.4)
            axes[0, 1].plot(t, a["min_bz_over_bias"], color=color, label=name)
            axes[1, 0].plot(t, a["ohmic_heat"], color=color, label=name)
            axes[1, 1].plot(t, a["peak_j_theta"]*1e-6, color=color, label=name)
        axes[0, 0].set(title="Center Bz · solid: conductor / dashed: vacuum", ylabel="Bz [mT]")
        axes[0, 1].set(title="Minimum Bz / Bbias within prescribed plasma", ylabel="Bz / Bbias")
        axes[1, 0].set(title="Accumulated Ohmic heat · no thermal feedback", ylabel="Heat [J]")
        axes[1, 1].set(title="Peak absolute induced toroidal current", ylabel="|Jθ| [MA/m²]")
        for ax in axes.flat:
            ax.set_xlabel("Time [µs]"); ax.grid(alpha=.18)
            ax.legend(framealpha=.1, labelcolor="#c9d5e5")
        for ax in axes[0]:
            ax.axhline(0, color="white", lw=.8, alpha=.4)
        peak = next(iter(results.values())).metadata["normalization"]["value_t"]*1e3
        fig.suptitle(f"PRESCRIBED WAVEFORM COMPARISON\nEqual peak |vacuum center perturbation| = {peak:g} mT · Stationary conductivity model", fontsize=16, weight="bold")
        fig.savefig(path, dpi=160)
        plt.close(fig)


def phase_diagram(rows, path):
    categories = ["no reversal", "transient reversal", "closed flux", "persistent closed flux (sampled)"]
    cmap = ListedColormap(["#334155", "#528cbf", "#e6ad69", "#75d8ba"])
    sigma = sorted({x["conductivity_s_m"] for x in rows})
    ratios = sorted({x["peak_external_over_bias"] for x in rows})
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(12, 4.8), layout="constrained")
        for ax, name in zip(axes, COLORS):
            matrix = np.zeros((len(sigma), len(ratios)))
            for row in rows:
                if row["waveform"] == name:
                    matrix[sigma.index(row["conductivity_s_m"]), ratios.index(row["peak_external_over_bias"])] = categories.index(row["summary"]["classification"])
            im = ax.imshow(matrix, origin="lower", cmap=cmap, norm=BoundaryNorm(np.arange(5)-.5, 4), aspect="auto")
            ax.set(xticks=range(len(ratios)), xticklabels=ratios, yticks=range(len(sigma)),
                   yticklabels=[f"{s/1e3:g}" for s in sigma], xlabel="Peak |external B| / |bias|", title=name.upper())
            for i in range(len(sigma)):
                for j in range(len(ratios)):
                    ax.text(j, i, str(int(matrix[i, j])), ha="center", va="center", color="white", weight="bold")
        axes[0].set_ylabel("Conductivity [kS/m]")
        fig.suptitle("EXPLORATORY PHASE MAP · fixed geometry and prescribed conductivity\nClosed flux ≠ equilibrium, stability or particle confinement", fontsize=14)
        cb = fig.colorbar(im, ax=axes, orientation="horizontal", fraction=.12, pad=.14, ticks=range(4))
        cb.ax.set_xticklabels(["0 · no reversal", "1 · reversal only", "2 · closed flux", "3 · persistent closed flux"], fontsize=8)
        fig.savefig(path, dpi=170)
        plt.close(fig)
