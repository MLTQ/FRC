"""Physical-reference material comparisons with explicit boundary and support errors."""
from pathlib import Path
import hashlib
import json
import numpy as np
from .surface_measurements import protocol,provenance
from .surface_sampling import fields,nodal
from .material_overlap import intersections,quadrature as overlap_quadrature
from .surface_time_control import case as extra_time_case


def read_study(directory,include_extra=True):
    directory=Path(directory);source=provenance();rows={};arrays={}
    refs=json.loads((directory/"references/metrics.json").read_text())
    if not refs["completed"] or refs["provenance"]!=source:
        raise ValueError("Incomplete or changed reference study")
    for row in refs["cases"]:
        path=directory/"references"/row["fields"]["file"]
        if not row["completed"] or row["provenance"]!=source or hashlib.sha256(path.read_bytes()).hexdigest()!=row["fields"]["sha256"]:
            raise ValueError("Changed or incomplete independent reference")
    for case in protocol()+([extra_time_case()] if include_extra else []):
        row=json.loads((directory/"dynamic"/f"{case.name}.json").read_text());path=directory/"dynamic"/row["fields"]["file"]
        if not row["completed"] or row["provenance"]!=source or row["configuration"]!=case.configuration():
            raise ValueError(f"Incomplete or mismatched surface case {case.name}")
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row["fields"]["sha256"]:
            raise ValueError("Changed field archive")
        with np.load(path,allow_pickle=False) as data:
            arrays[case.name]={k:data[k] for k in data.files}
        rows[case.name]=row
    return refs,rows,arrays


def compare(test,reference,arrays,quadrature=4):
    a,b=arrays[test["name"]],arrays[reference["name"]];c,r=test["configuration"],reference["configuration"]
    if c["kind"]!=r["kind"] or c["diffusivity"]!=r["diffusivity"] or c["angular_cells"]>r["angular_cells"]:
        raise ValueError("Comparisons need the same physics and a containing reference polygon")
    times=test["fields"]["times"][1:]
    if any(t not in reference["fields"]["times"] for t in times):raise ValueError("Mismatched output time")
    totals=[{} for _ in times];volume=0.
    def accumulate(target,key,w,x,y):
        if x.ndim==2:w=w[:,None]
        value=np.array([np.sum(w*(x-y)**2),np.sum(w*y*y)])
        target[key]=target.get(key,np.zeros(2))+value
    polygons=intersections(a,b)
    for points,area_weight in overlap_quadrature(polygons,quadrature):
        mass=2*np.pi*points[:,0]*area_weight;volume+=float(mass.sum()) # Initial rho=1.
        sa0=fields(a,c,"t0",points);sb0=fields(b,r,"t0",points)
        for saved,total in zip(times,totals):
            p=saved["prefix"];sa=fields(a,c,p,points);sb=fields(b,r,p,points)
            weight=mass*sa["volume_ratio"]
            accumulate(total,"velocity_relative_mass_L2",mass,sa["velocity"],sb["velocity"])
            accumulate(total,"total_current_relative_volume_L2",weight,sa["current"],sb["current"])
            accumulate(total,"local_heat_density_relative_volume_L2",weight,sa["heat_density"],sb["heat_density"])
            for key in ("magnetic","current","pressure","density"):
                accumulate(total,key+"_change_relative_volume_L2",weight,sa[key]-sa0[key],sb[key]-sb0[key])
    volume_error=volume/float(a["t0_volume"].sum())-1
    if abs(volume_error)>2e-10:raise ValueError("Overlap integration did not cover the complete test polygon")
    boundary=a["t0_reference_coordinates"][a["t0_boundary"]]
    normal=boundary/np.linalg.norm(boundary,axis=1)[:,None]
    samples=[]
    for saved,total in zip(times,totals):
        p=saved["prefix"]
        va=nodal(a,c,p,boundary,"velocity");vb=nodal(b,r,p,boundary,"velocity")
        row={"time_s":saved["time_s"],"velocity_rms_difference_m_s":float(np.sqrt(total["velocity_relative_mass_L2"][0]/volume)),
             "maximum_boundary_velocity_difference_m_s":float(np.linalg.norm(va-vb,axis=1).max()),
             "maximum_boundary_normal_velocity_difference_m_s":float(np.max(abs(np.sum((va-vb)*normal,axis=1))))}
        row.update({key:float(np.sqrt(value[0]/max(value[1],np.finfo(float).tiny))) for key,value in total.items()})
        samples.append(row)
    return {"test":test["name"],"reference":reference["name"],"quadrature":quadrature,
            "comparison_support_volume_relative_error":volume_error,
            "omitted_reference_initial_volume_fraction":1-float(a["t0_volume"].sum()/b["t0_volume"].sum()),"samples":samples}


def comparison_pairs():
    pairs=[(f"source_g{g:g}_{n}",f"source_g{g:g}_{2*n}") for g in (1.,1.5,2.) for n in (8,16)]
    pairs += [(f"frc_g{g:g}_16",f"frc_g{g:g}_32") for g in (1.5,2.)]
    pairs += [("source_g1_32","source_g1.5_32"),("source_g1.5_32","source_g2_32")]
    pairs += [("source_g2_16",f"source_g2_16_{s}") for s in ("angular","time","tolerance","matching","quadrature")]
    pairs += [("source_g2_32",f"source_g2_32_{s}") for s in ("angular","cfl")]
    pairs += [("source_g2_32_cfl","source_g2_32_cfl_quarter")]
    return pairs


def comparisons(rows,arrays):
    return [compare(rows[a],rows[b],arrays) for a,b in comparison_pairs()]
