"""Independent parent-species, local energy, and cylindrical kinematics checks."""
import numpy as np
from scipy.constants import elementary_charge as E,mu_0 as MU0
from plasma_bench.finite_mass_moments import SpeciesMasses,moments,reduced_ohm
from examples.moving_species_reference import MovingSpeciesReference,parent_sources,advection,divergence


def fixture():
    ref=MovingSpeciesReference();r,z=np.meshgrid(np.linspace(.001,.034,9),np.linspace(-.11,.11,11),indexing='ij')
    return ref,r,z,.23/ref.frequency


def test_exact_species_moments_and_moving_equation_reduction():
    ref,r,z,t=fixture();a=ref.fields(r,z,t);dr,dz,dt=ref.derivatives(r,z,t)
    force,heat,continuity=parent_sources(ref,r,z,t);mass=SpeciesMasses(ref.ion_mass,ref.electron_mass)
    m=moments(mass,a['n'],a['v'],a['j'],a['b'],a['pe'],a['pi'],a['e'],gamma_e=ref.gamma_e,gamma_i=ref.gamma_i)
    np.testing.assert_allclose(m['total_energy_flux'],ref.species_flux(r,z,t),rtol=3e-14,atol=1e-9)
    for key in ('i','e'):
        assert np.linalg.norm(continuity[key])<1e-12*np.linalg.norm(dt['n'])
    kinetic=.5*a['n']*(ref.ion_mass*np.sum(a['vi']**2,axis=-1)+ref.electron_mass*np.sum(a['ve']**2,axis=-1))
    np.testing.assert_allclose(m['bulk_kinetic']+m['current_kinetic'],kinetic,rtol=1e-15)
    grad=lambda name:np.stack((dr[name],np.zeros_like(r),dz[name]),axis=-1)
    conv=lambda x,y:advection(a[x],a[y],dr[y],dz[y],r)
    bulk=mass.total*a['n'][...,None]*(dt['v']+conv('v','v'))-np.cross(a['j'],a['b'])+grad('pe')+grad('pi')
    bulk+=mass.reduced*a['n'][...,None]*conv('w','w')
    np.testing.assert_allclose(bulk,force['i']+force['e'],rtol=3e-13,atol=3e-11)
    ohm=reduced_ohm(mass,a['n'],a['v'],a['w'],a['b'],grad('pe'),grad('pi'),dt['w'],conv('v','w'),conv('w','v'),conv('w','w'))
    source=(mass.beta*force['e']-mass.alpha*force['i'])/(a['n'][...,None]*E)
    np.testing.assert_allclose(ohm+source,a['e'],rtol=4e-13,atol=4e-11)
    np.testing.assert_allclose(dt['w'],dt['j']/(E*a['n'][...,None])-a['w']*(dt['n']/a['n'])[...,None],rtol=2e-14,atol=1e-4)


def test_independent_local_energy_flux_and_faraday():
    ref,r,z,t=fixture();a=ref.fields(r,z,t);dr,dz,dt=ref.derivatives(r,z,t)
    forces,heat,_=parent_sources(ref,r,z,t)
    rate=np.sum(a['b']*dt['b'],axis=-1)/MU0+dt['pe']/(ref.gamma_e-1)+dt['pi']/(ref.gamma_i-1)
    work=heat['i']+heat['e']
    for key,mass in (('i',ref.ion_mass),('e',ref.electron_mass)):
        rate+=.5*mass*dt['n']*np.sum(a['v'+key]**2,axis=-1)+mass*a['n']*np.sum(a['v'+key]*dt['v'+key],axis=-1)
        work+=np.sum(a['v'+key]*forces[key],axis=-1)
    def derivative(axis,h):
        values=[ref.species_flux(r+i*h if axis==0 else r,z+i*h if axis==1 else z,t) for i in (-2,-1,1,2)]
        return (values[0]-8*values[1]+8*values[2]-values[3])/(12*h)
    errors=[]
    for h in (1e-3*ref.radius,5e-4*ref.radius):
        div=divergence(ref.species_flux(r,z,t),derivative(0,h),derivative(1,h),r)
        errors.append(np.linalg.norm(rate+div-work)/np.linalg.norm(work))
    assert errors[1]<errors[0]/10 and errors[1]<1e-9,errors
    minus_curl=np.stack((dz['e'][...,1],dr['e'][...,2]-dz['e'][...,0],-dr['e'][...,1]-a['e'][...,1]/r),axis=-1)
    np.testing.assert_allclose(minus_curl,dt['b'],rtol=1e-13,atol=1e-10)
