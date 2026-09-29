"""Inspect frozen complete material-plasma archives without altering evolved state."""
from pathlib import Path
import hashlib
import json
import numpy as np
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.moving_advection import MaterialFlux
from .moving_cases import provenance
from .moving_measurements import protocol


def read_study(directory):
    directory=Path(directory);source=provenance();cases={};arrays={}
    static=json.loads((directory/"static.json").read_text())
    if not static["completed"] or static["provenance"]!=source:
        raise ValueError("Static source mismatch or incomplete study")
    for spec in protocol():
        row=json.loads((directory/"cases"/f"{spec.name}.json").read_text())
        if not row["completed"] or row["provenance"]!=source or row["configuration"]!=spec.configuration():
            raise ValueError(f"Incomplete or mismatched case: {spec.name}")
        path=directory/"cases"/row["fields"]["file"]
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row["fields"]["sha256"]:
            raise ValueError(f"Changed field archive: {spec.name}")
        with np.load(path,allow_pickle=False) as archive:
            arrays[spec.name]={key:archive[key] for key in archive.files}
        cases[spec.name]=row
    return static,cases,arrays


def domain(row):
    c=row["configuration"]
    return MovingDomain(c["radial_cells"],c["angular_cells"],radius=c["radius_m"],
                        matching_radius=c["matching_radius"],vacuum_layers=c["vacuum_layers"],
                        modes=c["modes"],quadrature=c["quadrature"])


def interior_flux_error(row,arrays):
    """Actual interpolated material-point psi; nodal invariance alone is insufficient."""
    d=domain(row);cells=d.cells[:d.fluid_cell_count];shape=d.shape
    x0=arrays["t0_coordinates"];a0=arrays["t0_total_material_a"]
    r0=np.einsum("qn,cn->cq",shape,x0[cells,0])
    psi0=r0*r0*np.einsum("qn,cn->cq",shape,a0[cells])
    initial=d.elements(x0)[0]
    weights=initial.weights[:len(cells)]*initial.points[:len(cells),:,0]
    norm=float(np.sum(weights*psi0**2));result=[]
    for saved in row["fields"]["times"]:
        p=saved["prefix"];x=arrays[p+"_coordinates"];a=arrays[p+"_total_material_a"]
        r=np.einsum("qn,cn->cq",shape,x[cells,0]);psi=r*r*np.einsum("qn,cn->cq",shape,a[cells])
        result.append({"time_s":saved["time_s"],"maximum_interpolated_psi_drift":float(np.max(abs(psi-psi0))),
                       "material_point_psi_relative_L2":float(np.sqrt(np.sum(weights*(psi-psi0)**2)/norm)) if norm else 0.})
    return result


def compare(test,reference,arrays):
    c,r=test["configuration"],reference["configuration"]
    ratio=r["radial_cells"]//c["radial_cells"]
    if ratio not in (1,2) or ratio*c["radial_cells"]!=r["radial_cells"]:
        raise ValueError("Comparison requires identical or nested doubled material meshes")
    d=domain(reference)
    select=np.r_[0,d.ids[ratio:d.nr+1:ratio,::ratio].ravel()]
    a,b=arrays[test["name"]],arrays[reference["name"]];rows=[]
    for saved in test["fields"]["times"][1:]:
        p=saved["prefix"]
        if saved not in reference["fields"]["times"]:
            raise ValueError("Compare identical physical saved times")
        v,w=a[p+"_velocity"],b[p+"_velocity"][select];mass=a[p+"_node_mass"]
        dv=np.sum(mass[:,None]*(v-w)**2)
        rows.append({"time_s":saved["time_s"],"velocity_mass_rms_difference_m_s":float(np.sqrt(dv/mass.sum())),
                     "velocity_relative_mass_L2":float(np.sqrt(dv/np.sum(mass[:,None]*w*w))),
                     "maximum_position_difference_m":float(np.linalg.norm(a[p+"_coordinates"]-b[p+"_coordinates"][select],axis=-1).max())})
    return {"test":test["name"],"reference":reference["name"],"samples":rows}


def expansion_pressure_error(row,arrays):
    d=domain(row);batch=d.elements(arrays["t10_coordinates"])[0];n=d.fluid_cell_count
    scale=np.sqrt(1+2*row["final"]["time_s"]**2)
    exact=(1-np.sum(batch.points[:n]**2,axis=-1)/scale**2)/scale**5
    weights=2*np.pi*batch.weights[:n]*batch.points[:n,:,0]
    p=(weights*exact).sum(axis=1)/weights.sum(axis=1)
    volume=arrays["t10_volume"]
    return float(np.sqrt(np.sum(volume*(p-arrays["t10_pressure"])**2)/np.sum(volume*p*p)))


def analyze(cases,arrays):
    comparisons=[compare(cases[f"source_{n}"],cases[f"source_{2*n}"],arrays) for n in (4,8)]
    comparisons += [compare(cases[name],cases["frc_16"],arrays) for name in
                    ("frc_far_matching","frc_modes4","frc_quadrature7")]
    comparisons += [compare(cases[f"{kind}_time_{dt}"],cases[f"{kind}_time_{dt/2}"],arrays)
                    for kind in ("expansion","frc","source") for dt in (.002,.001,.0005)]
    return {"comparisons":comparisons,
            "interpolated_material_flux":{name:interior_flux_error(row,arrays[name]) for name,row in cases.items()},
            "expansion_pressure_cell_average_error":{name:expansion_pressure_error(row,arrays[name])
                for name,row in cases.items() if row["configuration"]["kind"]=="expansion"}}


def nonuniform_advection():
    rows=[]
    for n in (4,8,16,32):
        d=MovingDomain(n,2*n);x0=d.fluid_reference;x=x0.copy()
        x[:,0]*=1+.2*x0[:,1];x[:,1]+=.1*x0[:,0]**2
        a0=.5*(1-.2*x0[:,0]**2+.1*x0[:,1])
        a=MaterialFlux(d,a0).coefficients(x);cells=d.cells[:d.fluid_cell_count]
        r0=d.shape@x0[cells,0].T;r=d.shape@x[cells,0].T
        before=r0*r0*(d.shape@a0[cells].T);after=r*r*(d.shape@a[cells].T)
        initial=d.elements(x0)[0];weights=(initial.weights[:len(cells)]*initial.points[:len(cells),:,0]).T
        rows.append({"radial_cells":n,"angular_cells":2*n,
                     "material_point_psi_relative_L2":float(np.sqrt(np.sum(weights*(after-before)**2)/np.sum(weights*before**2)))})
    return {"initial_a":"0.5 (1-0.2 r^2+0.1 z) T", "radial_map":"r_new=r (1+0.2 z)",
            "axial_map":"z_new=z+0.1 r^2", "radius_m":1., "quadrature":5,"rows":rows,
            "interpretation":"Prescribed nonuniform material deformation, not dynamical evolution; drift is numerical, not resistive"}
