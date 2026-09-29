"""Diagnose saved variational forces and a separately recovered volume J cross B."""
import numpy as np
from scipy.constants import mu_0 as MU0
from .surface_measurements import protocol
from .spherical_frc_reference import SphericalFRCReference


def forces(name,arrays,prefix="t10",time=.04):
    case=next(c for c in protocol() if c.name==name);s=case.build();a=arrays[name]
    s.state.coordinates=a[prefix+"_coordinates"].copy();s.state.velocity=a[prefix+"_velocity"].copy()
    s.flux.labels=a[prefix+"_flux_labels"].copy();s.pressure_invariant=a[prefix+"_pressure_invariant"].copy()
    s.cell_mass=a[prefix+"_cell_mass"].copy();s.node_mass=a[prefix+"_node_mass"].copy();s.time=time
    e=s.evaluate();d=s.domain;nc=d.fluid_cell_count;batch=e["magnetic"].elements[0]
    archived=a[prefix+"_induced_a"]
    mismatch=float(np.max(abs(e["magnetic"].induced-archived)))
    if mismatch>1e-10*max(float(np.max(abs(archived))),np.finfo(float).tiny):
        raise ValueError("Restored force state does not reproduce archived magnetic solution")
    pressure=np.zeros_like(s.state.coordinates)
    np.add.at(pressure,d.cells[:nc].ravel(),(e["pressure"][:,None,None]*d.volume_gradients(e["magnetic"].elements)).reshape(-1,2))
    pressure[d.axis,0]=0.
    magnetic=e["force"]-pressure
    br,bz=s.field.fields(batch,e["magnetic"].induced,e["bias"])
    current=a[prefix+"_current_quadrature"];weight=2*np.pi*batch.weights[:nc]*batch.points[:nc,:,0]
    local=np.einsum("cq,cqn,cqa->cna",weight,batch.geometry_shape[:nc],
                    current[...,None]*np.stack((bz[:nc],-br[:nc]),axis=-1))
    direct=np.zeros_like(pressure);np.add.at(direct,d.cells[:nc].ravel(),local.reshape(-1,2));direct[d.axis,0]=0.
    ids=d.ids[:d.nr+1,d.nt//2];mass=s.node_mass[ids]
    result={"reference_radius_m":d.fluid_reference[ids,0],"radial_velocity_m_s":s.state.velocity[ids,0],
            "pressure_radial_acceleration_m_s2":pressure[ids,0]/mass,
            "magnetic_radial_acceleration_m_s2":magnetic[ids,0]/mass,
            "weak_current_radial_acceleration_m_s2":direct[ids,0]/mass,
            "net_radial_acceleration_m_s2":e["force"][ids,0]/mass,
            "restored_induced_max_difference_t":mismatch}
    return {key:value.tolist() if isinstance(value,np.ndarray) else value for key,value in result.items()}


def initial_compatibility():
    """Check the smooth p=0 trace against the fixture's initial Joule source."""
    theta=np.linspace(.2,np.pi-.2,17);epsilon=1e-5;step=1e-6;rows=[]
    for name,bias in (("frc",np.sqrt(MU0)),("source",np.sqrt(MU0)*np.exp(-.5))):
        ref=SphericalFRCReference(bias=bias)
        r=(1-epsilon)*np.sin(theta);z=(1-epsilon)*np.cos(theta)
        br_plus=ref.fields(r,z+step)[1];br_minus=ref.fields(r,z-step)[1]
        bz_plus=ref.fields(r+step,z)[2];bz_minus=ref.fields(r-step,z)[2]
        recovered=((br_plus-br_minus)-(bz_plus-bz_minus))/(2*step*MU0*(1-epsilon))
        exact=-7.5*bias/MU0*np.sin(theta)
        error=float(np.max(abs(recovered-exact))/np.max(abs(exact)))
        if error>1e-7:raise ValueError("Independent boundary-current derivative check failed")
        rate=(2/3)*(.1*MU0)*exact**2
        rows.append({"fixture":name,"radius_m":1.,"gamma":5/3,"physical_diffusivity_m2_s":.1,
                     "initial_bias_t":float(bias),"theta_radians":theta.tolist(),"material_current_trace_a_m2":exact.tolist(),
                     "initial_pressure_source_pa_s":rate.tolist(),"equatorial_pressure_source_pa_s":float(rate[8]),
                     "one_sided_maxwell_current_relative_error":error,
                     "interpretation":"Conditional incompatibility with smooth material p=0, bounded div(v), continuous B and no other heat/stress terms; not a general weak-solution nonexistence claim"})
    return rows
