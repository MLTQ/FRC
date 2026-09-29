"""Computed resistive current/pressure fields and independent accuracy figures."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from .moving_plot import field_panel
from .moving_analysis import domain


def error_figure(directory,decay,rows):
    fig,axes=plt.subplots(2,3,figsize=(14,8.4),layout="constrained")
    mesh=[r for r in decay["cases"] if r["name"].startswith("mesh_")]
    n=[r["configuration"]["radial_cells"] for r in mesh]
    for key,label in (("field_relative_volume_L2","B: Q1 gradient"),("flux_relative_volume_L2","poloidal flux"),
                      ("current_relative_volume_L2","current density")):
        axes[0,0].loglog(n,[r["history"][-1][key] for r in mesh],"o-",label=label)
    axes[0,0].set(title="Independent sphere decay, t=0.1 s",ylabel="Relative volume L2 error")
    for key,label in (("local_heat_relative_volume_L2","local Ohmic heat"),
                      ("maximum_boundary_current_scaled_error","surface current / reference"),
                      ("axis_a_relative_error","axis field")):
        axes[0,1].loglog(n,[abs(r["history"][-1][key]) for r in mesh],"o-",label=label)
    axes[0,1].set(title="Heating, interface and regular axis",ylabel="Relative / scaled analytic error")
    stiff=decay["stiff_control"]
    axes[0,2].loglog([r["tolerance"] for r in stiff["adaptive"]],
                     [r["amplitude_absolute_error"] for r in stiff["adaptive"]],"o-")
    axes[0,2].text(.04,.96,"Raw single step: amplitude −1/3\nExact decay: exp(−4) > 0",transform=axes[0,2].transAxes,
                   va="top",fontsize=9,bbox=dict(facecolor="white",edgecolor=".8",alpha=.9))
    axes[0,2].set(title="Adaptive control resolves false reversal",xlabel="Local diffusion tolerance",ylabel="Amplitude error / initial amplitude")
    cells=[4,8,16,32]
    axes[1,0].plot(cells,[rows[f"frc_{i}"]["final"]["maximum_speed_m_s"] for i in cells],"o-")
    axes[1,0].set(title="Initialized FRC: resistive relaxation",xlabel="Radial cells",ylabel="Final maximum speed [m/s]",xticks=cells)
    axes[1,1].plot(cells,[rows[f"source_{i}"]["final"]["maximum_speed_m_s"] for i in cells],"o-",label="D=0.1 m²/s")
    axes[1,1].plot([16,32],[rows[f"source_wide_{i}"]["final"]["maximum_speed_m_s"] for i in (16,32)],"s-",label="D=1 m²/s")
    axes[1,1].set(title="Driven flow: different physical D",xlabel="Radial cells",ylabel="Final maximum speed [m/s]",xticks=cells)
    times=[rows[name] for name in ("source_time_0.001","source_8","source_time_0.00025")]
    axes[1,2].loglog([r["configuration"]["max_timestep"] for r in times],
                     [abs(r["final"]["response_energy_residual_j"]) for r in times],"o-")
    axes[1,2].set(title="Moving source: whole-response error",xlabel="Outer time cap [s]",ylabel="|Response-energy residual| [J]")
    for ax in axes.flat:
        ax.grid(alpha=.25,which="both")
    for ax in axes[0,:2]:
        ax.set(xlabel="Radial cells",xticks=n);ax.set_xticklabels(n)
        ax.xaxis.set_minor_formatter(NullFormatter());ax.legend(fontsize=8)
    axes[1,1].legend(fontsize=9)
    axes[1,2].set_xticks([.00025,.0005,.001],labels=["0.00025","0.0005","0.001"])
    axes[1,2].xaxis.set_minor_formatter(NullFormatter())
    fig.suptitle("Resistive moving plasma / infinite vacuum — numerical verification",fontsize=15)
    fig.savefig(Path(directory)/"resistive_verification.png",dpi=170);plt.close(fig)


def field_figure(directory,rows,arrays):
    fig,axes=plt.subplots(2,2,figsize=(12,8.2),layout="constrained")
    row=rows["source_32"];a=arrays["source_32"];prefix="t10"
    pressure=field_panel(axes[0,0],row,a,prefix,"Computed plasma pressure")
    fig.colorbar(pressure,ax=axes[0,0],label="p [Pa]",shrink=.8)
    current=field_panel(axes[0,1],row,a,prefix,"Computed toroidal current density")
    d=domain(row);e=d.elements(a[prefix+"_coordinates"])[0];nc=d.fluid_cell_count
    weight=e.weights[:nc]*e.points[:nc,:,0]
    average=np.sum(weight*a[prefix+"_current_quadrature"],axis=1)/weight.sum(axis=1)
    current.set_array(average);current.set_cmap("RdBu_r");bound=np.max(abs(average));current.set_clim(-bound,bound)
    fig.colorbar(current,ax=axes[0,1],label="Jθ [A/m²]",shrink=.8)
    for n in (8,16,32):
        r=rows[f"source_{n}"];h=r["history"];ar=arrays[r["name"]]
        axes[1,0].plot([v["time_s"] for v in h],[v["ohmic_heat_j"] for v in h],label=f"nr={n}")
        ids=ar["t10_boundary"];initial=ar["t0_coordinates"][ids];velocity=ar["t10_velocity"][ids]
        angle=np.arctan2(initial[:,0],initial[:,1])*180/np.pi
        outward=np.sum(velocity*initial,axis=1)/np.linalg.norm(initial,axis=1)
        axes[1,1].plot(angle,outward,label=f"nr={n}")
    axes[1,0].set(title="Positive physical Ohmic heating",xlabel="Time [s]",ylabel="Cumulative heat [J]")
    axes[1,1].set(title="Surface speed along initial radius",xlabel="Initial polar angle [degrees]",ylabel="Outward velocity component [m/s]")
    for ax in axes[1]:
        ax.grid(alpha=.25);ax.legend(fontsize=9)
    fig.suptitle("Initialized FRC under a prescribed bias ramp, D=0.1 m²/s, t=0.04 s\nComputed fields and flux contours; no formation or stability claim",fontsize=13)
    fig.savefig(Path(directory)/"resistive_fields.png",dpi=170);plt.close(fig)
