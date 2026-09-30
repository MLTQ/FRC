"""Independent smooth forces and original-current induction fixture."""
import numpy as np
from scipy.constants import electron_mass
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.regular_field_map import RegularFieldMap
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from .variable_current_reference import VariableCurrentReference
from .finite_mass_manufactured import cell_average


def force_regular(r,z,ref):
    y=(r/ref.radius)**2;phase=np.pi*z/ref.half_length
    return np.stack(np.broadcast_arrays(np.exp(y)*(1+.2*np.sin(phase)),
        np.cos(y)*(1+.3*np.cos(phase)),np.exp(.2*y)*(1+.1*np.sin(phase))),axis=-1)


def force(r,z,ref):
    value=force_regular(r,z,ref);value[...,:2]*=r[...,None]
    return value


def initialize(spec):
    ref=VariableCurrentReference(contrast=spec.get('contrast',.8),number_density=spec.get('density',1e18))
    ne=spec['elements'];degree=spec.get('degree',5);points=spec.get('points',4*ne+8)
    d=RegularCoordinateFluid(ref.radius,2*ref.half_length,ne,points,ref.ion_mass+electron_mass,degree=degree)
    r,z=np.broadcast_arrays(d.r,d.z);f=ref.fields(r,z);n=f['n'];y=(r/ref.radius)**2;phase=np.pi*z/ref.half_length
    vr_regular=ref.number_density*ref.speed/ref.radius*(1-y)**2*np.cos(phase)/n
    omega=.1*ref.speed/ref.radius*(1+.2*y)*np.cos(phase)
    state=np.stack((n,n*vr_regular,n*omega,f['particle'][...,2],n*np.log(2.),n*np.log(1.5)),axis=-1);state[0,0,:,1:3]=0
    nr=spec.get('field_nr',max(8,2*ne));nz=spec.get('field_nz',nr)
    mesh=MHDMesh(ref.radius,ref.half_length,nr,nz)
    geo=CompatibleGeometry(mesh,degree=spec.get('field_degree',5),quadrature_order=spec.get('quadrature',12))
    mapping=RegularFieldMap(d,geo);model=RegularCurrent(mapping,n,ref.ion_mass)
    psi=ref.fields(*mesh.corners)['psi'];psi[[0,-1]]=0.;psi[:,-1]=psi[:,0]
    g=cell_average(mesh,lambda r,z:ref.fields(r,z)['g'],order=12)
    field=HarmonicMagnetic(mesh,psi,g,ref.axial_bias)
    return ref,d,state,model,field
