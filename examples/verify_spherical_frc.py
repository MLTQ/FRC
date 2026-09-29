"""Render a prescribed analytic free-boundary FRC reference, never formation."""
from dataclasses import asdict
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import mu_0 as MU0
from .spherical_frc_reference import SphericalFRCReference
from plasma_bench.geometry import Grid
from plasma_bench.topology import analyze_flux


def quadrature(reference,order=24):
    a,b=reference.radius,reference.bias
    nodes,w=leggauss(order);mu=nodes[None,:]
    s=.5*a*(nodes[:,None]+1)
    _,br,bz,p,_=reference.fields(s*np.sqrt(1-mu*mu),s*mu)
    volume=2*np.pi*s*s*(.5*a*w[:,None])*w[None,:]
    inside=np.sum((br*br+(bz-b)**2)/(2*MU0)*volume)
    total_inside=np.sum((br*br+bz*bz)/(2*MU0)*volume)
    pressure=np.sum(p*volume)
    x=.5*(nodes[:,None]+1);s=a/x
    _,br,bz,_,_=reference.fields(s*np.sqrt(1-mu*mu),s*mu)
    volume=2*np.pi*a**3/x**4*(.5*w[:,None])*w[None,:]
    outside=np.sum((br*br+(bz-b)**2)/(2*MU0)*volume)
    actual={"pressure_integral_j":float(pressure),"inside_total_magnetic_energy_j":float(total_inside),
            "induced_field_energy_j":float(inside+outside),"exterior_induced_field_energy_j":float(outside)}
    exact=reference.integrals()
    return {"gauss_order":order,"values":actual,
            "relative_errors":{key:abs(value/exact[key]-1) for key,value in actual.items()}}


def verify(directory="outputs/free_boundary_reference"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from plasma_bench.visualization import STYLE
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    ref=SphericalFRCReference()
    r=np.linspace(0,1.35*ref.radius,301)
    z=np.linspace(-1.5*ref.radius,1.5*ref.radius,601)
    rr,zz=np.meshgrid(r,z,indexing="ij")
    psi,br,bz,p,j=ref.fields(rr,zz)
    data=(bz/ref.bias,p/(ref.bias**2/MU0),j*MU0*ref.radius/ref.bias)
    titles=("Axial field and closed flux","Pressure · exterior exactly zero","Toroidal current · vacuum exactly zero")
    labels=("Bz / imposed bias","p / (B0² / mu0)","Jθ / [B0 / (mu0 a)]")
    with plt.rc_context(STYLE):
        fig,axes=plt.subplots(1,3,figsize=(16,3.8),layout="constrained")
        theta=np.linspace(0,np.pi,401)
        for ax,value,title,label,cmap in zip(axes,data,titles,labels,("coolwarm","viridis","magma_r")):
            im=ax.pcolormesh(z/ref.radius,r/ref.radius,value,shading="auto",cmap=cmap)
            fig.colorbar(im,ax=ax,label=label,shrink=.8)
            levels=np.linspace(-3*ref.bias*ref.radius**2/16,0,9)[1:-1]
            ax.contour(zz/ref.radius,rr/ref.radius,psi,levels=levels,colors="#f6f7fb",linewidths=.6,alpha=.65)
            ax.plot(np.cos(theta),np.sin(theta),color="#f6f7fb",linewidth=1.6,label="Spherical separatrix")
            ax.plot(0,1/np.sqrt(2),"o",color="#ffcc75",markersize=5,label="Magnetic O-ring")
            ax.set(title=title,xlabel="z / a",ylabel="r / a",aspect="equal",xlim=(-1.5,1.5),ylim=(0,1.35))
        axes[0].legend(fontsize=7,loc="upper right",labelcolor="#c9d5e5")
        fig.suptitle("ANALYTIC SPHERICAL FRC / MAGNETIC VACUUM\nPrescribed static equilibrium · no time evolution or formation result",fontsize=15,weight="bold")
        fig.savefig(directory/"spherical_frc.png",dpi=150);plt.close(fig)
    topology=[]
    for n in (24,48,96,192):
        grid=Grid(1.4*ref.radius,1.5*ref.radius,n+1,2*n+1)
        diagnostic=analyze_flux(grid,ref.fields(*grid.mesh)[0])
        topology.append({"grid":asdict(grid),"diagnostics":diagnostic,
                         "separatrix_volume_relative_error":abs(diagnostic["separatrix_volume_m3"]/ref.integrals()["separatrix_volume_m3"]-1)})
    files=[Path(__file__),Path("examples/spherical_frc_reference.py")]
    files += [Path("plasma_bench")/f"{name}.py" for name in ("geometry","numerics","topology")]
    record={"interpretation":"Prescribed analytic equilibrium; not formation, stability, confinement or a transport closure",
            "parameters":asdict(ref),"exact_integrals":ref.integrals(),
            "quadrature":[quadrature(ref,n) for n in (8,16,24)],
            "topology_refinement":topology,"python":platform.python_version(),
            "versions":{name:importlib.metadata.version(name) for name in ("numpy","scipy","contourpy","matplotlib")},
            "source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    (directory/"metrics.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    print(json.dumps(record["quadrature"],indent=2),flush=True)
    return record


if __name__=="__main__":
    verify()
