"""Shared smooth conservative mass, momentum and enthalpy face quadrature."""
from functools import lru_cache
import numpy as np
from scipy.constants import mu_0 as MU0
from .annular_polynomial import radial_matrix,axial_samples
from .annular_cell_polynomial import radial_cell_matrix,axial_cell_samples
from .mhd_state import NonphysicalState

POWERS=(0,1,0,2,0,0,0,1,0,0)


@lru_cache(maxsize=32)
def matrices(mesh,order):
    nodes,weights=np.polynomial.legendre.leggauss(order)
    faces={power:radial_matrix(mesh,mesh.rf,power=power) for power in (0,1,2)}
    cells=[{power:radial_cell_matrix(mesh,.5*(node+1),power=power) for power in (0,1,2)} for node in nodes]
    for matrix in [*faces.values(),*(m for row in cells for m in row.values())]:matrix.setflags(write=False)
    return nodes,weights,faces,cells


def variables(state,evaluation,parameters):
    u=state.conserved;q=evaluation.primitive
    pe=(parameters.gamma_e-1)*evaluation.components['electron']
    pi=(parameters.gamma_i-1)*evaluation.components['ion']
    return np.stack((u[...,0],u[...,1],u[...,2],u[...,4],pe,pi,
                     evaluation.components['current'],q[...,4],q[...,5],state.g),axis=-1)


def radial(values,matrices):
    return np.stack([matrices[power]@values[...,index] for index,power in enumerate(POWERS)],axis=-1)


def physical(value,radius,axis,parameters):
    rho,mr,mz,angular,pe,pi,kj,br,bz,g=np.moveaxis(value,-1,0)
    if not np.isfinite(value).all() or np.any(rho<=0) or np.any(pe<=0) or np.any(pi<=0):
        raise NonphysicalState('Invalid smooth reconstructed density/species pressure; no floor or fallback')
    vr,vz=mr/rho,mz/rho
    vt=np.divide(angular,rho*radius,out=np.zeros_like(angular),where=radius!=0)
    bt=radius*g;magnetic=(br*br+bz*bz+bt*bt)/(2*MU0)
    mass=mr if axis==0 else mz;bn=br if axis==0 else bz
    enthalpy=parameters.gamma_e/(parameters.gamma_e-1)*pe+parameters.gamma_i/(parameters.gamma_i-1)*pi
    specific=.5*(vr*vr+vz*vz+vt*vt)+(kj+enthalpy)/rho
    flux=np.stack((mass,vr*mass-br*bn/MU0,vz*mass-bz*bn/MU0,
                   specific*mass,angular/rho*mass-radius*bt*bn/MU0),axis=-1)
    flux[...,axis+1]+=pe+pi+magnetic
    hoop=pe+pi+magnetic+rho*vt*vt-bt*bt/MU0
    return flux,hoop


def bulk_fluxes(mesh,state,evaluation,parameters,*,quadrature=5):
    """No artificial diffusion: bounded smooth central reference, not a shock flux."""
    if type(quadrature) is not int or quadrature<5:raise ValueError('At least fifth-order face quadrature required')
    nodes,weights,face_matrices,cell_matrices=matrices(mesh,quadrature)
    values=variables(state,evaluation,parameters)
    # Normal momentum is reconstructed as a conserved quantity, not rho times
    # independently limited velocity. Transverse integration retains its mean.
    radial_average=radial(values,face_matrices);radial_average[[0,-1],:,1]=0.
    axial_average=axial_samples(values,0.)
    br,bz=mesh.face_fields(state.psi);radial_average[...,7]=br
    axial_average[...,8]=bz[:,:-1]
    radial_average[-1,:,9]=0.  # Explicit Btheta(R)=0 boundary traction.
    fr=np.zeros((mesh.nr+1,mesh.nz,5));fz=np.zeros((mesh.nr,mesh.nz,5));source=np.zeros((mesh.nr,mesh.nz))
    for node,weight,cell_matrix in zip(nodes,weights,cell_matrices):
        rr=mesh.r[:,None]+.5*mesh.dr*node
        rf_values=axial_cell_samples(radial_average,.5*(node+1))
        radial_flux,_=physical(rf_values,mesh.rf[:,None],0,parameters)
        axial_flux,_=physical(radial(axial_average,cell_matrix),rr,1,parameters)
        fr+=.5*weight*radial_flux
        fz+=(.5*weight*rr*mesh.dr/mesh.radial_measure[:,None])[...,None]*axial_flux
        rv=radial(values,cell_matrix)
        for z,wz in zip(nodes,weights):
            _,hoop=physical(axial_cell_samples(rv,.5*(z+1)),rr,0,parameters)
            source+=.25*weight*wz*mesh.dr/mesh.radial_measure[:,None]*hoop
    fr*=2*np.pi*mesh.rf[:,None,None]*mesh.dz
    fz*=2*np.pi*mesh.radial_measure[:,None,None]
    fz=np.concatenate((fz,fz[:,:1]),axis=1)
    return fr,fz,source
