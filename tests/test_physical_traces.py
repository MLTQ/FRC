"""Smooth physical polynomials and exact manufactured one-sided jumps."""
import numpy as np
from examples.strong_wall_reference import initialize
from examples.physical_trace_sampling import TraceSampler


def setup():
    _,d,u,s,b=initialize(dict(kind='contact',elements=2,points=8,field_nr=12,field_nz=8))
    return d,u,s.geometry.mesh


def test_exact_polynomial_current_and_magnetic_face_values():
    d,u,m=setup();r=m.rf[2:-1];z=m.zf[:-1]
    x=.5*m.rf**2;theta=(m.nr-1)*m.nz
    jr=np.broadcast_to(2+3*(x[:-1]+x[1:])[:,None]/2,(m.nr,m.nz))
    jz=np.broadcast_to((4+5*x+6*x*x)[:,None],(m.nr+1,m.nz))
    current=np.r_[jr.ravel(),np.zeros(theta),jz.ravel()]
    psi=np.broadcast_to((x+2*x*x)[:,None],(m.nr+1,m.nz));magnetic=np.r_[psi.ravel(),np.full(m.nr*m.nz,.7)]
    for shift in (0,-1):
        sample=TraceSampler(d,m,r,z,np.arange(2,m.nr)+shift,np.arange(m.nz)+shift)
        j,gradient=sample.electric(current,True);b=sample.magnetic(magnetic)
        assert np.allclose(j[...,0],sample.R*(2+1.5*sample.R**2),rtol=2e-12,atol=0)
        assert np.allclose(j[...,2],4+2.5*sample.R**2+1.5*sample.R**4,rtol=2e-12,atol=0)
        assert np.allclose(np.trace(gradient,axis1=-2,axis2=-1),4+6*sample.R**2,rtol=2e-10,atol=1e-10)
        assert np.max(abs(b[...,0]))<1e-13
        assert np.allclose(b[...,1],.7*sample.R,rtol=2e-12,atol=0)
        assert np.allclose(b[...,2],1+2*sample.R**2,rtol=2e-12,atol=0)


def test_degree_six_flux_has_known_distinct_one_sided_degree_five_derivatives():
    d,u,m=setup();face=6;x=.5*m.rf**2;scale=x[-1]
    psi=np.broadcast_to((x/scale)[:,None]**6,(m.nr+1,m.nz));magnetic=np.r_[psi.ravel(),np.zeros(m.nr*m.nz)]
    actual=[]
    for cell in (face-1,face):
        sample=TraceSampler(d,m,m.rf[[face]],np.array([m.z[0]]),np.array([cell]),np.array([0]))
        value=sample.magnetic(magnetic)[0,0,2];actual.append(value)
        ids=np.arange(cell-2,cell+4);others=x[ids[ids!=face]]
        expected=6*x[face]**5/scale**6-np.prod(x[face]-others)/scale**6
        assert np.isclose(value,expected,rtol=2e-13)
    assert not np.isclose(*actual,rtol=1e-6)


def test_radial_fluid_interface_retains_both_density_traces():
    d,u,m=setup();u[...,0]=2e20;u[1,...,0]=3e20;face=m.nr//2
    for cell,expected in [(face-1,2e20),(face,3e20)]:
        sample=TraceSampler(d,m,m.rf[[face]],m.z[:2],np.array([cell]),np.arange(2))
        assert np.allclose(sample.scalar(u[...,0]),expected,rtol=2e-14)
