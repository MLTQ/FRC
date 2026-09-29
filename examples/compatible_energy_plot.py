"""Scientific comparison of compatible energy, source budgets and RK defects."""
import numpy as np


def plot(metrics,directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    from plasma_bench.visualization import STYLE
    rows = {r["name"]:r for r in metrics["cases"]}
    old = {r["name"]:r for r in metrics["historical_local"]["cases"]}
    colors = {"Historical cell / Riemann / RK2":"#ffcc75", "Staggered / variational / RK3":"#72ddc4"}
    with plt.rc_context(STYLE):
        fig,axes = plt.subplots(2,3,figsize=(16,9.5),layout="constrained")
        for label,data,prefix in ((next(iter(colors)),old,"shock_local_"),
                                  (list(colors)[1],rows,"shock_")):
            shock = [data[f"{prefix}{n}"] for n in (24,48,96,192,384)
                     if data.get(f"{prefix}{n}",{}).get("completed",False)]
            x = [r["radial_cells"] for r in shock]
            axes[0,0].loglog(x,[r["errors"]["Bz_L1_over_upstream_B"] for r in shock],"o-",color=colors[label],label=label)
            axes[0,1].loglog(x,[abs(r["diagnostics"]["response_energy_balance_residual_j"]) for r in shock],
                             "o-",color=colors[label],label=label)
            pulse = [data[f"pulse_{n}"] for n in (16,32,64) if data.get(f"pulse_{n}",{}).get("completed",False)]
            axes[1,0].loglog([r["radial_cells"] for r in pulse],
                             [100*abs(r["diagnostics"]["response_residual_over_source_work"]) for r in pulse],
                             "o-",color=colors[label],label=label)
        axes[0,0].set(title="Original shock · interior magnetic error",xlabel="Radial cells",ylabel="Bz L1 / upstream field")
        axes[0,1].set(title="Independent full-space energy error",xlabel="Radial cells",ylabel="Absolute response residual [J]")
        axes[1,0].set(title="Same coil pulse · source budget",xlabel="Radial cells",ylabel="Absolute response error / applied work [%]")
        probes=[r for r in metrics.get("boundary_order_probe",[]) if r["completed"] and r["cfl"]==.2]
        probes.sort(key=lambda r:r["radial_cells"])
        if probes:
            x=[r["radial_cells"] for r in probes]
            axes[0,0].loglog(x,[r["errors"]["Bz_L1_over_upstream_B"] for r in probes],"x--",color="#e69ce9",label="Experimental boundary strip / RK3")
            axes[0,1].loglog(x,[abs(r["diagnostics"]["response_energy_balance_residual_j"]) for r in probes],"x--",color="#e69ce9",label="Experimental boundary strip / RK3")
        for order,color in ((2,"#ffcc75"),(3,"#72ddc4")):
            names = ("shock_rk2_192","shock_rk2_half_cfl") if order==2 else ("shock_192","shock_half_cfl")
            chosen = [rows[n] for n in names if rows.get(n,{}).get("completed",False)]
            axes[0,2].plot([r["parameters"]["cfl"] for r in chosen],
                           [abs(r["diagnostics"]["response_energy_balance_residual_j"]) for r in chosen],
                           "o-",color=color,label=f"Compatible RK{order}")
        fine_probes=[r for r in metrics.get("boundary_order_probe",[]) if r["completed"] and r["radial_cells"]==192]
        if fine_probes:
            fine_probes.sort(key=lambda r:r["cfl"])
            axes[0,2].plot([r["cfl"] for r in fine_probes],[abs(r["diagnostics"]["response_energy_balance_residual_j"]) for r in fine_probes],
                           "x--",color="#e69ce9",label="Experimental boundary strip / RK3")
        axes[0,2].set(title="nr=192 · measured timestep sensitivity",xlabel="CFL",ylabel="Absolute response residual [J]",xticks=[.1,.2],yscale="log")
        for audit in metrics.get("time_audit",[]):
            if audit["radial_cells"] != 192:
                continue
            data = audit["snapshots"]
            axes[1,1].semilogy([r["time_s"] for r in data],[abs(r["measured_response_change_j"]) for r in data],
                               "o-",color="#aa9cfc",label="Measured energy error")
            axes[1,1].semilogy([r["time_s"] for r in data],[abs(r["accumulated_temporal_prediction_j"]) for r in data],
                               "x--",color="#ffcc75",label="Independent RK prediction")
        axes[1,1].set(title="Energy error explained by RK stages",xlabel="Time [s]",ylabel="Absolute cumulative error [J]")
        for kind,color in (("advection","#72ddc4"),("diffusion","#aa9cfc")):
            selected = [rows[f"{kind}_{n}"] for n in (12,24,48) if rows.get(f"{kind}_{n}",{}).get("completed",False)]
            axes[1,2].loglog([r["radial_cells"] for r in selected],[r["flux_relative_L2"] for r in selected],
                             "o-",color=color,label=kind.capitalize())
        axes[1,2].set(title="Independent smooth-solution checks",xlabel="Radial cells",ylabel="Flux relative L2 error")
        for ax in axes.flat:
            ax.grid(alpha=.2)
            handles,labels = ax.get_legend_handles_labels()
            if handles and ax is not axes[0,1]:
                ax.legend(fontsize=7,labelcolor="#c9d5e5")
        for ax in (axes[0,0],axes[0,1]):
            ax.set_xticks([24,48,96,192,384],labels=["24","48","96","192","384"])
            ax.xaxis.set_minor_formatter(NullFormatter())
        axes[1,0].set_xticks([16,32,64],labels=["16","32","64"])
        axes[1,2].set_xticks([12,24,48],labels=["12","24","48"])
        axes[1,0].xaxis.set_minor_formatter(NullFormatter());axes[1,2].xaxis.set_minor_formatter(NullFormatter())
        fig.suptitle("COMPATIBLE MAGNETIC ENERGY / VACUUM WORK\nFixed-chamber verification · unchanged shock and sources · no energy repair",fontsize=15,weight="bold")
        fig.savefig(directory/"compatible_energy.png",dpi=150);plt.close(fig)
