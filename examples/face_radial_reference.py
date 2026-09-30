"""Independent physical initialization of the staggered radial momentum state."""
import numpy as np
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.face_radial_fluid import FaceRadialFluid
from .log_fluid_reference import LogFluidReference,average


def initial(driver,reference):
    m=driver.mesh
    density=average(m,lambda r,z:reference.fields(r,z)['n'])
    momentum=average(m,lambda r,z:reference.fields(r,z)['momentum'][...,1:])
    q,w=np.polynomial.legendre.leggauss(16)
    r=m.rf[1:-1,None,None];z=m.z[None,:,None]+.5*m.dz*q
    radial=m.rf[1:-1,None]*np.sum(.5*w*reference.fields(r,z)['momentum'][...,0],axis=-1)
    r,z=np.meshgrid(m.rf,m.zf[:-1],indexing='ij');f=reference.fields(r,z)
    return driver.initialize(density,radial,momentum,f['electron'],f['ion'])


def initialize(spec,reference=None):
    kind=spec.get('kind','contact')
    ref=reference or LogFluidReference(kind=kind,direction=spec.get('direction',1),gamma_i=5/3 if kind=='acoustic' else 1.4)
    nr=spec['nr'];g=CompatibleGeometry(MHDMesh(ref.radius,ref.half_length,nr,2*nr),degree=spec.get('degree',5),quadrature_order=spec.get('quadrature_order',12))
    solver=FaceRadialFluid(g,ref.ion_mass,gamma_e=ref.gamma_e,gamma_i=ref.gamma_i)
    return ref,solver,initial(solver,ref)
