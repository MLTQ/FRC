"""Computed graded-mesh profiles, force cancellation and independent accuracy."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def profiles(directory,rows,arrays,force_rows):
    fig,axes=plt.subplots(2,3,figsize=(17,10),layout="constrained")
    names=["source_g1_32","source_g1.5_32","source_g2_16","source_g2_32","source_g2_32_angular"]
    labels=["uniform 32 × 64","grade 1.5, 32 × 64","grade 2, 16 × 32","grade 2, 32 × 64","grade 2, 32 × 128"]
    for name,label in zip(names,labels):
        a=arrays[name];c=rows[name]["configuration"];nt=c["angular_cells"];nr=c["radial_cells"]
        ids=np.r_[0,1+np.arange(nr)*(nt+1)+nt//2]
        x=a["t0_reference_coordinates"][ids,0];p="t10"
        axes[0,0].plot(x,a[p+"_velocity"][ids,0],".-",label=label)
        cells=np.arange(nr)*nt+nt//2;centers=a["t0_reference_coordinates"][a["t0_cells"][cells]].mean(axis=1)[:,0]
        axes[0,1].plot(centers,a[p+"_pressure"][cells],".-",label=label)
        axes[0,2].plot(centers,a[p+"_density"][cells],".-",label=label)
        axes[1,0].plot(x,a[p+"_coordinates"][ids,0]*a[p+"_current_over_radius"][ids]/1000,".-",label=label)
        boundary=a[p+"_boundary"];initial=a["t0_reference_coordinates"][boundary]
        normal=initial/np.linalg.norm(initial,axis=1)[:,None]
        axes[1,1].plot(np.linspace(0.,180,nt+1),np.sum(a[p+"_velocity"][boundary]*normal,axis=1),label=label)
    for ax,title,ylabel in zip(axes.ravel()[:4],["Equatorial material flow","Equatorial cell pressure","Equatorial cell density","Recovered toroidal current"],
                               ["Radial velocity [m/s]","Pressure [Pa]","Density [kg/m³]","Jθ [kA/m²]"]):
        ax.set(title=title,ylabel=ylabel,xlabel="Initial material radius [m]",xlim=(.8,1.005));ax.grid(alpha=.25)
    axes[0,0].legend(fontsize=8)
    axes[1,1].set(title="Actual surface motion",xlabel="Initial polar angle [degrees]",ylabel="Outward velocity component [m/s]")
    axes[1,1].grid(alpha=.25)
    f=force_rows["source_g2_32"]
    for key,label in [("pressure","pressure"),("magnetic","magnetic virtual work"),("weak_current","recovered volume J × B"),("net","net")]:
        axes[1,2].plot(f["reference_radius_m"],f[key+"_radial_acceleration_m_s2"],".-",label=label)
    axes[1,2].set(title="Force cancellation: grade 2, 32 × 64",xlabel="Initial material radius [m]",ylabel="Radial acceleration [m/s²]",xlim=(.9,1.005))
    axes[1,2].legend(fontsize=8);axes[1,2].grid(alpha=.25)
    fig.suptitle("Initialized FRC under a bias ramp: unchanged D=0.1 m²/s, t=0.04 s\nComputed surface layers; no formation or convergence claim",fontsize=17)
    fig.savefig(Path(directory)/"surface_profiles.png",dpi=140);plt.close(fig)


def accuracy(directory,refs,comparisons):
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout="constrained")
    for g,color in [(1.5,"tab:blue"),(2.,"tab:orange")]:
        decay=[r for r in refs["cases"] if r["name"].startswith(f"decay_g{g:g}_") and not r["name"].endswith("angular")]
        n=[r["configuration"]["radial_cells"] for r in decay]
        axes[0,0].loglog(n,[r["history"][-1]["field_relative_volume_L2"] for r in decay],"o-",color=color,label=f"grade {g:g}")
        for key,style,label in [("current_relative_volume_L2","o-","current"),("local_heat_relative_volume_L2","s--","local heat")]:
            axes[0,1].loglog(n,[r["history"][-1][key] for r in decay],style,color=color,label=f"{label}, grade {g:g}")
        gas=[r for r in refs["cases"] if r["name"].startswith(f"expansion_g{g:g}_")]
        axes[1,0].loglog([r["configuration"]["radial_cells"] for r in gas],[r["errors"]["velocity_relative_mass_L2"] for r in gas],"o-",color=color,label=f"grade {g:g}")
    for g in (1.,1.5,2.):
        values=[next(r for r in comparisons if r["test"]==f"source_g{g:g}_{n}" and r["reference"]==f"source_g{g:g}_{2*n}") for n in (8,16)]
        axes[1,1].plot([8,16],[100*max(s["velocity_relative_mass_L2"] for s in r["samples"]) for r in values],"o-",label=f"grade {g:g}")
    titles=["Independent sphere decay: magnetic field","Independent sphere decay: current and heat","Independent free gas expansion","Driven flow: adjacent meshes"]
    ylabels=["Relative volume L2 error","Relative volume L2 error","Velocity relative mass L2 error","Velocity difference [% of finer reference]"]
    for i,(ax,title,ylabel) in enumerate(zip(axes.ravel(),titles,ylabels)):
        ax.set(title=title,xlabel="Radial cells" if i<3 else "Coarser radial cells",ylabel=ylabel)
        ax.set_xticks([4,8,16] if i==2 else [8,16] if i==3 else [8,16,32]);ax.set_xticklabels([4,8,16] if i==2 else [8,16] if i==3 else [8,16,32])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(alpha=.25,which="both");ax.legend(fontsize=9)
    fig.suptitle("Surface refinement: independent errors and actual moving-state differences",fontsize=15)
    fig.savefig(Path(directory)/"surface_accuracy.png",dpi=140);plt.close(fig)
