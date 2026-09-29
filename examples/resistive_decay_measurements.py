"""Continuum sphere-decay profiles and independent stiff-time accuracy controls."""
from pathlib import Path
import hashlib
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.linalg import eigh
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.moving_field import MovingField
from plasma_bench.material_resistivity import MaterialResistivity
from plasma_bench.resistive_substeps import integrate_diffusion
from .spherical_decay_reference import SphericalDecay
from .resistive_records import provenance
from .verification_records import write_record


def diagnostics(o,a,heat,ref,time):
    d=o.domain;nc=d.fluid_cell_count;e=o.elements[0]
    r,z=e.points[:nc,:,0],e.points[:nc,:,1];weight=2*np.pi*e.weights[:nc]*r
    u=o.induced(a);q=o.current_over_radius(a);j=o.current_quadrature(q)
    exact=ref.fields(r,z,time);br,bz=MovingField.fields(e,u)
    value=np.einsum("cn,cqn->cq",a[e.field_ids[:nc]],e.field_shape[:nc])
    target_heat=np.sum(weight*ref.heat_density(r,z,time),axis=1);volume=weight.sum(axis=1)
    berror=np.sqrt(np.sum(weight*((br[:nc]-exact["Br"])**2+(bz[:nc]-exact["Bz"])**2))/np.sum(weight*(exact["Br"]**2+exact["Bz"]**2)))
    jerror=np.sqrt(np.sum(weight*(j-exact["Jtheta"])**2)/np.sum(weight*exact["Jtheta"]**2))
    psierror=np.sqrt(np.sum(weight*(r*r*value-exact["psi"])**2)/np.sum(weight*exact["psi"]**2))
    heaterror=np.sqrt(np.sum((heat-target_heat)**2/volume)/np.sum(target_heat**2/volume)) if time else 0.
    boundary=d.fluid_reference[d.boundary];ab=ref.radial(ref.radius,time)[0]
    jboundary=boundary[:,0]*q[d.boundary];jtarget=ref.wavenumber**2/MU0*boundary[:,0]*ab
    boundary_scale=3*abs(ref.central_a)/(MU0*ref.radius)*np.exp(-ref.diffusivity*ref.wavenumber**2*time)
    energy=float(np.pi/MU0*np.dot(u,o.stiffness@u));expected=ref.integrals(time)
    return {"time_s":time,"field_relative_volume_L2":float(berror),"current_relative_volume_L2":float(jerror),
            "flux_relative_volume_L2":float(psierror),"local_heat_relative_volume_L2":float(heaterror),
            "maximum_boundary_current_scaled_error":float(np.max(abs(jboundary-jtarget))/boundary_scale),
            "axis_a_relative_error":float(a[0]/ref.fields(0.,0.,time)["a"]-1),
            "induced_energy_j":energy,"analytic_induced_energy_j":expected["induced_energy_j"],
            "energy_relative_error":energy/expected["induced_energy_j"]-1,
            "ohmic_heat_j":float(heat.sum()),"analytic_ohmic_heat_j":expected["ohmic_heat_j"]},u,q,j


def decay_case(name,n,dt=.002,matching=2.,modes=8,quadrature=5,radial_mode=1,tolerance=1e-4):
    ref=SphericalDecay(mode=radial_mode)
    d=MovingDomain(n,2*n,matching_radius=matching,vacuum_layers=round(n*(matching-1)),modes=modes,quadrature=quadrature)
    o=MaterialResistivity(d,d.fluid_reference,ref.diffusivity*MU0)
    a=ref.fields(*d.fluid_reference.T)["a"];heat=np.zeros(d.fluid_cell_count)
    row={"name":name,"completed":True,"configuration":{"radial_cells":n,"angular_cells":2*n,"radius_m":1.,
         "matching_radius_m":matching,"vacuum_layers":round(n*(matching-1)),"modes":modes,"quadrature":quadrature,
         "physical_diffusivity_m2_s":ref.diffusivity,"physical_resistivity_ohm_m":ref.diffusivity*MU0,
         "central_initial_a_t":ref.central_a,"radial_mode":radial_mode,"outer_time_cap_s":dt,
         "end_time_s":.1,"diffusion_tolerance":tolerance,"held_fixed":"all material geometry; induction component only"},
         "history":[],"accepted_substeps":0,"rejected_trials":0,"maximum_accepted_error":0.,"maximum_substep_balance_j":0.}
    time=0.;arrays={}
    for index,end in enumerate(np.linspace(0.,.1,11)):
        while end-time>1e-14:
            step=min(dt,end-time)
            result=integrate_diffusion(o,a,step,tolerance=tolerance)
            a=result.material_a;heat+=result.cell_heat;time+=step
            row["accepted_substeps"]+=result.accepted_substeps;row["rejected_trials"]+=result.rejected_trials
            row["maximum_accepted_error"]=max(row["maximum_accepted_error"],result.maximum_accepted_error)
            row["maximum_substep_balance_j"]=max(row["maximum_substep_balance_j"],result.maximum_balance_residual)
        diag,u,q,j=diagnostics(o,a,heat,ref,float(end));row["history"].append(diag)
        if index in (0,5,10):
            arrays.update({f"t{index}_material_a":a.copy(),f"t{index}_induced":u,f"t{index}_current_over_radius":q,
                           f"t{index}_current_quadrature":j,f"t{index}_cell_heat":heat.copy()})
    arrays.update(coordinates=d.fluid_reference,cells=d.cells[:d.fluid_cell_count],
                  quadrature_coordinates=o.elements[0].points[:d.fluid_cell_count],
                  quadrature_volume=2*np.pi*o.elements[0].weights[:d.fluid_cell_count]*o.elements[0].points[:d.fluid_cell_count,:,0])
    row["response_residual_j"]=row["history"][-1]["induced_energy_j"]+heat.sum()-row["history"][0]["induced_energy_j"]
    return row,arrays


def stiff_control():
    d=MovingDomain(4,8);o=MaterialResistivity(d,d.fluid_reference,.1*MU0)
    n=d.fluid_node_count;k=o.stiffness.toarray();m=o.mass.toarray()
    s=k[:n,:n]-k[:n,n:]@np.linalg.solve(k[n:,n:],k[n:,:n]);values,vectors=eigh(s,m)
    initial=vectors[:,0]*np.sqrt(MU0)/vectors[0,0];duration=4/(.1*values[0]);raw=o.step(initial,duration)
    rows=[]
    for tol in (1e-3,1e-4,1e-5,1e-6):
        r=integrate_diffusion(o,initial,duration,tolerance=tol)
        amp=float((initial@m@r.material_a)/(initial@m@initial))
        rows.append({"tolerance":tol,"amplitude":amp,"amplitude_absolute_error":abs(amp-np.exp(-4)),
                     "accepted_substeps":r.accepted_substeps,"rejected_trials":r.rejected_trials,
                     "maximum_accepted_error":r.maximum_accepted_error,"heat_j":float(r.cell_heat.sum())})
    return {"interpretation":"Exact independently diagonalized semidiscrete low mode; isolates time error, not spatial accuracy",
            "duration_s":float(duration),"exact_amplitude":float(np.exp(-4)),"raw_midpoint_amplitude":float((initial@m@raw.material_a)/(initial@m@initial)),
            "raw_balance_residual_j":raw.balance_residual,"adaptive":rows}


def run(directory="outputs/resistive_moving_verification/decay"):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);source=provenance();rows=[]
    cases=[(f"mesh_{n}",n,{}) for n in (4,8,16,32,64)]
    cases += [("time_half",32,{"dt":.001}),("tolerance",32,{"tolerance":1e-5}),
              ("far_matching",16,{"matching":3.}),("near_matching",16,{"matching":1.5}),
              ("modes4",16,{"modes":4}),("quadrature7",16,{"quadrature":7}),
              ("mode2_16",16,{"radial_mode":2}),("mode2_32",32,{"radial_mode":2})]
    for name,n,options in cases:
        row,arrays=decay_case(name,n,**options);path=directory/f"{name}.npz"
        np.savez_compressed(path,**arrays)
        row.update(provenance=source,fields={"file":path.name,"sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "times_s":[0.,.05,.1],"units":{"material_a":"T","induced":"T, modal coefficients last",
             "current_over_radius":"A/m^3","current_quadrature":"A/m^2","cell_heat":"J",
             "coordinates":"m, r/z","quadrature_coordinates":"m, r/z","quadrature_volume":"m^3"}})
        write_record(directory/f"{name}.json",row);rows.append(row)
        print("Finished decay",name,flush=True)
    if provenance()!=source:
        raise RuntimeError("Source changed during decay study")
    write_record(directory/"metrics.json",{"completed":True,"provenance":source,"cases":rows,"stiff_control":stiff_control()})


if __name__=="__main__":
    run()
