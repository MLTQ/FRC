"""Scientific error plots and computed material/FRC field visualizations."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.ticker import NullFormatter
from .moving_analysis import domain


def plot_errors(directory,static,cases):
    fig,axes=plt.subplots(2,3,figsize=(14,8.3),layout="constrained")
    for scale in (.9,1.,1.1):
        rows=[r for r in static["homologous"] if r["scale"]==scale]
        n=[r["radial_cells"] for r in rows]
        axes[0,0].loglog(n,[abs(r["energy_relative_error"]) for r in rows],"o-",label=f"scale {scale:g}")
        axes[0,1].loglog(n,[abs(r["denergy_dscale_relative_error"]) for r in rows],"o-",label=f"scale {scale:g}")
    axes[0,0].set(title="Magnetized sphere: energy",ylabel="Relative analytic error")
    axes[0,1].set(title="Magnetized sphere: virtual force",ylabel="Relative analytic error")
    n=[4,8,16,32]
    axes[0,2].loglog(n,[abs(r["spurious_axial_force_n"]) for r in static["translation"]],"o-",color="tab:red")
    axes[0,2].set(title="Translated sphere: spurious force",ylabel="|Net axial force| [N]")
    rows=[cases[f"expansion_{i}"] for i in n]
    axes[1,0].loglog(n,[r["final_errors"]["velocity_relative_mass_L2"] for r in rows],"o-",label="velocity relative L2")
    axes[1,0].loglog(n,[r["final_errors"]["maximum_boundary_position_error_m"] for r in rows],"s-",label="surface position / 1 m")
    axes[1,0].set(title="Exact gas expansion, t=0.1 s",ylabel="Error")
    rows=[cases[f"frc_{i}"] for i in n]
    axes[1,1].loglog(n,[r["final"]["maximum_speed_m_s"] for r in rows],"o-",label="maximum speed")
    axes[1,1].loglog(n,[r["initial_errors"]["initial_or_current_acceleration_mass_rms_m_s2"]*.04 for r in rows],"s--",label="initial RMS acceleration × 0.04 s")
    axes[1,1].set(title="Initialized FRC: artificial motion",ylabel="Speed [m/s]")
    rows=[cases[f"source_time_{dt}"] for dt in (.001,.0005,.00025)]
    axes[1,2].loglog([r["configuration"]["max_timestep"] for r in rows],
                     [abs(r["final"]["response_energy_residual_j"]) for r in rows],"o-")
    axes[1,2].set(title="Changing bias: source-work error",xlabel="Time cap [s]",ylabel="|Response residual| [J]")
    for ax in axes.flat:
        ax.grid(alpha=.25,which="both")
    for ax in list(axes[0])+list(axes[1,:2]):
        ax.set(xlabel="Radial cells",xticks=n);ax.set_xticklabels(n)
        ax.xaxis.set_minor_formatter(NullFormatter())
    axes[1,2].set_xticks([.00025,.0005,.001],labels=["0.00025","0.0005","0.001"])
    axes[1,2].xaxis.set_minor_formatter(NullFormatter())
    for ax in (axes[0,0],axes[0,1],axes[1,0],axes[1,1]):
        ax.legend(fontsize=8)
    fig.suptitle("Moving plasma / infinite vacuum — independent numerical verification",fontsize=15)
    fig.savefig(Path(directory)/"convergence.png",dpi=170);plt.close(fig)


def field_panel(ax,row,arrays,prefix,title):
    d=domain(row);cells=arrays[prefix+"_cells"];x=arrays[prefix+"_coordinates"]
    a=arrays[prefix+"_total_material_a"];pressure=arrays[prefix+"_pressure"]
    polygons=x[cells][:,:,[1,0]]
    collection=PolyCollection(polygons,array=pressure,cmap="viridis",edgecolors="none")
    ax.add_collection(collection)
    # Subsample the actual Q1 interpolant, including r², rather than smoothing nodal psi.
    v=np.linspace(0.,1.,6);xi,tau=np.meshgrid(v,v,indexing="ij")
    shape=np.stack(((1-xi)*(1-tau),xi*(1-tau),xi*tau,(1-xi)*tau),axis=-1)
    points=np.einsum("ijn,cna->cija",shape,x[cells]);value=points[...,0]**2*np.einsum("ijn,cn->cij",shape,a[cells])
    count=points.shape[0];base=np.arange(count)*36;triangles=[]
    for i in range(5):
        for j in range(5):
            lo=base+6*i+j
            triangles.extend((np.column_stack((lo,lo+6,lo+7)),np.column_stack((lo,lo+7,lo+1))))
    triangles=np.concatenate(triangles);points=points.reshape(-1,2)
    edge1=points[triangles[:,1]]-points[triangles[:,0]];edge2=points[triangles[:,2]]-points[triangles[:,0]]
    signed=edge1[:,0]*edge2[:,1]-edge1[:,1]*edge2[:,0]
    triangles=triangles[abs(signed)>1e-18]
    low=float(value.min());levels=np.linspace(low*.95,low*.05,10)
    ax.tricontour(points[:,1],points[:,0],triangles,value.ravel(),levels=levels,colors="white",linewidths=.6)
    boundary=x[arrays[prefix+"_boundary"]]
    ax.plot(boundary[:,1],boundary[:,0],color="tab:red",lw=1.2,label="actual material surface")
    ax.set(xlim=(-1.04,1.04),ylim=(0,1.04),xlabel="z [m]",ylabel="r [m]",title=title,aspect="equal")
    ax.set_facecolor("#e8edf2");ax.text(.02,.96,"vacuum",transform=ax.transAxes,va="top",fontsize=9)
    return collection


def plot_fields(directory,cases,arrays):
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout="constrained")
    for ax,name,prefix,title in ((axes[0,0],"frc_32","t0","Initialized FRC, t=0"),
                               (axes[0,1],"source_16","t10","Prescribed bias ramp, t=0.04 s")):
        collection=field_panel(ax,cases[name],arrays[name],prefix,title)
        fig.colorbar(collection,ax=ax,shrink=.8,label="Plasma pressure [Pa]")
    for name in ("frc_8","frc_16","frc_32"):
        h=cases[name]["history"]
        axes[1,0].plot([r["time_s"] for r in h],[r["maximum_speed_m_s"] for r in h],label=name.replace("frc_","nr="))
    axes[1,0].set(title="Initialized equilibrium: numerical drift",xlabel="Time [s]",ylabel="Maximum artificial speed [m/s]")
    for n in (4,8,16,32):
        r=cases[f"expansion_{n}"];a=arrays[r["name"]];x=a["t10_coordinates"];x0=a["t0_coordinates"];idx=a["t10_boundary"]
        theta=np.arctan2(x0[idx,0],x0[idx,1])
        error=np.linalg.norm(x[idx],axis=-1)-np.sqrt(1.02)
        axes[1,1].plot(theta*180/np.pi,error*1e3,label=f"nr={n}")
    axes[1,1].set(title="Gas expansion: surface-radius error at 0.1 s",xlabel="Polar angle [degrees]",ylabel="Radial error [mm]")
    for ax in axes[1]:
        ax.legend(fontsize=9);ax.grid(alpha=.25)
    fig.suptitle("Computed material motion and flux contours — initialized references, not formation",fontsize=14)
    fig.savefig(Path(directory)/"material_fields.png",dpi=170);plt.close(fig)
