"""Figure: simulated NUCTE-III separatrix evolution against the published case-1 values."""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from examples.nucte_compare import summarize, reference

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]  # ordinal blue, step 250 -> 700
REFERENCE = "#eb6834"  # categorical slot 2 marks the published value


def style(ax, title, ylabel, xlabel="time (µs)"):
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    ax.set_xlabel(xlabel, color=MUTED, fontsize=9)
    ax.set_ylabel(ylabel, color=MUTED, fontsize=9)
    ax.grid(color=GRID, linewidth=.8)
    ax.tick_params(colors=MUTED, labelsize=8)
    for spine in ax.spines.values():
        spine.set_visible(False)


def main(runs, out):
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.2), constrained_layout=True)
    labels = []
    for k, run in enumerate(runs):
        setup, rf, zf, rows = summarize(run)
        if k == 0:
            case_label, target = reference(setup)
        label = Path(run).name
        labels.append(label)
        t = np.array([r["time_s"] for r in rows])*1e6
        get = lambda key, scale=1.: np.array([np.nan if r.get(key) is None else r[key]*scale for r in rows])
        marker = "o" if k == 0 else "s"
        if k == 0:
            ax = axes[0, 0]
            picks = [i for i, ti in enumerate(t) if round(ti) in (2, 4, 6, 8, 10, 15)]
            for colour, i in zip(RAMP, picks):
                prof = rows[i].get("profile")
                if prof is None:
                    continue
                prof = np.array([np.nan if x is None or not np.isfinite(x) else x for x in prof], dtype=float)
                ax.plot(zf*100, prof*100, color=colour, linewidth=2, label=f"{t[i]:.0f} µs")
            ax.axhline(rf[-1]*100, color=MUTED, linewidth=1, linestyle=":")
            ax.text(zf[0]*100, rf[-1]*100+.3, "tube wall", color=MUTED, fontsize=8)
            style(ax, f"Separatrix r_s(z), {label}", "r_s (cm)", "z (cm)")
            ax.legend(frameon=False, fontsize=8, labelcolor=INK, ncol=2)
        series = [(axes[0, 1], get("separatrix_mid_m", 100), target["separatrix_mid_m"]*100, "Midplane separatrix radius", "cm"),
                  (axes[0, 2], get("length_m", 100), target["separatrix_length_m"]*100, "Reversed-region length", "cm"),
                  (axes[1, 0], get("density_m3"), target["density_m3"], "Mean density inside separatrix", "m⁻³"),
                  (axes[1, 1], get("total_temperature_ev"), target["total_temperature_ev"], "Te+Ti inside separatrix", "eV"),
                  (axes[1, 2], get("trapped_flux_wb", 1e3), None, "Trapped poloidal flux", "mWb")]
        for ax, y, ref, title, unit in series:
            ax.plot(t, y, color=RAMP[2+2*k], linewidth=2, marker=marker, markersize=4, label=label)
            if k == 0 and ref is not None:
                ax.axhline(ref, color=REFERENCE, linewidth=2, linestyle="--", label=f"NUCTE-III {case_label} (equilibrium)")
            style(ax, title, unit)
    for ax in (axes[0, 1], axes[0, 2], axes[1, 0], axes[1, 1], axes[1, 2]):
        ax.legend(frameon=False, fontsize=8, labelcolor=INK)
    fig.suptitle(f"NUCTE-III formation: coarse single-fluid simulation vs published {case_label} (Takahashi 2000)",
                 x=.01, ha="left", color=INK, fontsize=12)
    fig.savefig(out, dpi=140, facecolor="#fcfcfb")
    print(out)


if __name__ == "__main__":
    main(sys.argv[1:-1], sys.argv[-1])
