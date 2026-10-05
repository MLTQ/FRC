"""Manufactured original local energy transport with a closed cylindrical wall."""
from types import SimpleNamespace
import numpy as np
from examples.strong_wall_reference import initialize
from examples.physical_cell_budget import measure
from plasma_bench.regular_current import RegularCurrent


def test_original_polynomial_compression_has_local_particle_entropy_energy_balance():
    _,d,u,s,b=initialize(dict(kind='contact',elements=2,points=8,field_nr=12,field_nz=8))
    n=d.thermal.n0;a=1700.;r=d.r;radius=d.radius;f=1-r*r/radius**2
    u[...,0]=n;u[...,1]=n*a*f;u[...,2:4]=0.
    u[...,4:]=n*np.log([2.,.6]);rate=np.zeros_like(u)
    divv=2*a*(1-2*r*r/radius**2);nt=-n*divv
    rate[...,0]=nt;rate[...,1]=nt*a*f-n*a*a*f*(1-3*r*r/radius**2)
    rate[...,4:]=nt[...,None]*np.log([2.,.6])
    m=s.geometry.mesh
    psi_t=-a*b.axial_bias*m.corners[0]**2*(1-m.corners[0]**2/radius**2)
    psi_t=np.broadcast_to(psi_t,(m.nr+1,m.nz+1)).copy();psi_t[[0,-1]]=0.
    bt=s.geometry.pack(psi_t,np.zeros((m.nr,m.nz)))
    et=a*b.axial_bias*m.rf[1:-1]*(1-m.rf[1:-1]**2/radius**2)
    electric=np.r_[np.zeros(m.nr*m.nz),np.broadcast_to(et[:,None],(m.nr-1,m.nz)).ravel(),np.zeros((m.nr+1)*m.nz)]
    model=RegularCurrent(s.mapping,u[...,0],s.ion_mass)
    manufactured=SimpleNamespace(fluid=rate,magnetic=bt,current=model.current(bt),electric=electric)
    result,arrays=measure(d,s,u,b,manufactured,12)
    for name in ('particles','electron_entropy','ion_entropy','energy'):
        assert result['channels'][name]['cell_balance_relative_L2']<2e-10
    assert result['channels']['current']['cell_balance_absolute_L2']==0.
