"""Hash-checked material-coordinate comparisons of completed resistive records."""
from pathlib import Path
import hashlib
import json
import numpy as np
from .resistive_records import provenance
from .resistive_dynamic_measurements import protocol


def read_study(directory):
    directory=Path(directory);source=provenance();rows={};arrays={}
    decay=json.loads((directory/"decay/metrics.json").read_text())
    if not decay["completed"] or decay["provenance"]!=source:
        raise ValueError("Incomplete or mismatched decay study")
    for r in decay["cases"]:
        path=directory/"decay"/r["fields"]["file"]
        if not r["completed"] or r["provenance"]!=source or hashlib.sha256(path.read_bytes()).hexdigest()!=r["fields"]["sha256"]:
            raise ValueError("Changed or incomplete decay case")
    for case in protocol():
        r=json.loads((directory/"dynamic"/f"{case.name}.json").read_text())
        path=directory/"dynamic"/r["fields"]["file"]
        if not r["completed"] or r["provenance"]!=source or r["configuration"]!=case.configuration():
            raise ValueError(f"Incomplete or mismatched dynamic case {case.name}")
        if hashlib.sha256(path.read_bytes()).hexdigest()!=r["fields"]["sha256"]:
            raise ValueError(f"Changed dynamic field archive {case.name}")
        with np.load(path,allow_pickle=False) as archive:
            arrays[case.name]={key:archive[key] for key in archive.files}
        rows[case.name]=r
    return decay,rows,arrays


def material_samples(arrays,prefix,coarse_n,fine_n,quadrature=5,subdivisions=1):
    """Evaluate actual Q1 B and recovered J at corresponding reference-sector points."""
    ratio=fine_n//coarse_n
    if ratio not in (1,2) or coarse_n*ratio!=fine_n:
        raise ValueError("Only matching or nested doubled material sectors are supported")
    v,w=np.polynomial.legendre.leggauss(quadrature)
    nodes=np.concatenate([(i+(v+1)/2)/subdivisions for i in range(subdivisions)])
    weights=np.tile(w/(2*subdivisions),subdivisions)
    xi,tau=np.meshgrid(nodes,nodes,indexing="ij");xi=xi.ravel();tau=tau.ravel()
    radial,angular=np.divmod(np.arange(coarse_n*2*coarse_n),2*coarse_n)
    ir=np.floor(ratio*xi).astype(int);it=np.floor(ratio*tau).astype(int)
    cells=(ratio*radial[:,None]+ir[None,:])*(2*fine_n)+ratio*angular[:,None]+it[None,:]
    xi=ratio*xi-ir;tau=ratio*tau-it
    shape=np.stack(((1-xi)*(1-tau),xi*(1-tau),xi*tau,(1-xi)*tau),axis=-1)
    deriv=np.stack((np.stack((-(1-tau),1-tau,tau,-tau),axis=-1),
                    np.stack((-(1-xi),-xi,xi,1-xi),axis=-1)),axis=-1)
    ids=arrays[prefix+"_cells"][cells];coordinates=arrays[prefix+"_coordinates"][ids]
    x=np.einsum("qn,cqna->cqa",shape,coordinates)
    jac=np.einsum("qnb,cqna->cqab",deriv,coordinates)
    a=arrays[prefix+"_total_material_a"][ids]
    gradient=np.einsum("cqnb,cqba->cqna",np.broadcast_to(deriv,(len(cells),*deriv.shape)),np.linalg.inv(jac))
    value=np.einsum("qn,cqn->cq",shape,a);da=np.einsum("cqna,cqn->cqa",gradient,a)
    r=x[...,0];magnetic=np.stack((-r*da[...,1],2*value+r*da[...,0]),axis=-1)
    current=r*np.einsum("qn,cqn->cq",shape,arrays[prefix+"_current_over_radius"][ids])
    weight=(weights[:,None]*weights[None,:]).ravel()[None,:]*2*np.pi*r*np.linalg.det(jac)*ratio**2
    return magnetic,current,weight


def restrict_cells(arrays,prefix,coarse_n,fine_n,key):
    ratio=fine_n//coarse_n;volume=arrays[prefix+"_volume"]
    v=volume.reshape(coarse_n,ratio,2*coarse_n,ratio).sum(axis=(1,3)).ravel()
    value=arrays[prefix+"_"+key]
    if key=="cell_ohmic_heat":
        integrated=value
    else:
        integrated=value*volume
    return integrated.reshape(coarse_n,ratio,2*coarse_n,ratio).sum(axis=(1,3)).ravel()/v


def compare(test,reference,arrays):
    c,r=test["configuration"],reference["configuration"]
    n,nf=c["radial_cells"],r["radial_cells"];ratio=nf//n
    if ratio not in (1,2) or n*ratio!=nf or c["diffusivity"]!=r["diffusivity"]:
        raise ValueError("Accuracy comparison must retain physical D on matching/nested meshes")
    a,b=arrays[test["name"]],arrays[reference["name"]]
    id_grid=np.arange(1,1+nf*(2*nf+1)).reshape(nf,2*nf+1)
    select=np.r_[0,id_grid[ratio-1::ratio,::ratio].ravel()]
    b0,j0,_=material_samples(a,"t0",n,n,subdivisions=ratio)
    br0,jr0,_=material_samples(b,"t0",n,nf,subdivisions=ratio)
    result=[]
    def relative(weight,value,reference):
        return float(np.sqrt(np.sum(weight*(value-reference)**2)/max(np.sum(weight*reference**2),np.finfo(float).tiny)))
    for saved in test["fields"]["times"][1:]:
        if saved not in reference["fields"]["times"]:
            raise ValueError("Comparisons require identical saved physical times")
        p=saved["prefix"];mass=a[p+"_node_mass"];volume=a[p+"_volume"]
        velocity=a[p+"_velocity"];vr=b[p+"_velocity"][select]
        bt,jt,weight=material_samples(a,p,n,n,subdivisions=ratio)
        br,jr,_=material_samples(b,p,n,nf,subdivisions=ratio)
        row={"time_s":saved["time_s"],"velocity_relative_mass_L2":relative(mass[:,None],velocity,vr),
             "velocity_rms_difference_m_s":float(np.sqrt(np.sum(mass[:,None]*(velocity-vr)**2)/mass.sum())),
             "maximum_position_difference_m":float(np.linalg.norm(a[p+"_coordinates"]-b[p+"_coordinates"][select],axis=-1).max()),
             "magnetic_change_relative_volume_L2":relative(weight[...,None],bt-b0,br-br0),
             "current_relative_volume_L2":relative(weight,jt,jr),
             "current_change_relative_volume_L2":relative(weight,jt-j0,jr-jr0)}
        for key in ("density","pressure"):
            reference_now=restrict_cells(b,p,n,nf,key);reference_initial=restrict_cells(b,"t0",n,nf,key)
            row[key+"_change_relative_volume_L2"]=relative(volume,a[p+"_"+key]-a["t0_"+key],reference_now-reference_initial)
        row["local_heat_density_relative_volume_L2"]=relative(volume,a[p+"_cell_ohmic_heat"]/volume,
                                                               restrict_cells(b,p,n,nf,"cell_ohmic_heat"))
        result.append(row)
    return {"test":test["name"],"reference":reference["name"],"samples":result}


def comparisons(rows,arrays):
    pairs=[(f"{kind}_{n}",f"{kind}_{2*n}") for kind in ("frc","source") for n in (4,8,16)]
    pairs += [("source_wide_16","source_wide_32")]
    for kind in ("frc","source"):
        pairs += [(f"{kind}_time_0.001",f"{kind}_8"),(f"{kind}_8",f"{kind}_time_0.00025"),
                  (f"{kind}_far_matching",f"{kind}_16"),(f"{kind}_tolerance",f"{kind}_16")]
    pairs += [("source_modes4","source_16"),("source_quadrature7","source_16")]
    return [compare(rows[a],rows[b],arrays) for a,b in pairs]
