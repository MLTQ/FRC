"""Independent polynomial, continuity, compression and physical entropy checks."""
import numpy as np
import pytest
from scipy.constants import atomic_mass,elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.variable_compatible_current import VariableCompatibleCurrent
from plasma_bench.compatible_log_thermal import CompatibleLogThermal
from plasma_bench.compatible_bulk_thermal import CompatibleBulkThermal
from plasma_bench.compatible_particle_flux import CompatibleParticleFlux,ParticleFaces
from plasma_bench.particle_flux_basis import primitive_basis,difference
from examples.finite_mass_manufactured import cell_average
from examples.compatible_field_norms import direct_basis
from examples.variable_current_reference import VariableCurrentReference


def initialize(nr=8,gamma=5/3):
    ref=VariableCurrentReference();mesh=MHDMesh(ref.radius,ref.half_length,nr,2*nr)
    n=cell_average(mesh,lambda r,z:ref.fields(r,z)['n'],order=12)
    model=VariableCompatibleCurrent(mesh,n,ref.ion_mass,quadrature_order=12)
    thermal=CompatibleLogThermal(model,gamma=gamma)
    return ref,model,thermal,CompatibleParticleFlux(model.density_map)


def test_flux_primitive_independent_polynomials_and_face_traces():
    edges=np.linspace(0,1,13)**2;points=np.concatenate([np.linspace(a,b,9)[:-1] for a,b in zip(edges[:-1],edges[1:])])
    polynomial=np.polynomial.Polynomial([.4,-.8,.6,.3,-.2,.1,.07])
    R,H=primitive_basis(edges,points);D=difference(edges)
    faces=polynomial(edges);derivative=polynomial.deriv()
    averages=np.diff(faces)/np.diff(edges)
    np.testing.assert_allclose(H@faces,polynomial(points),rtol=3e-12,atol=2e-13)
    np.testing.assert_allclose(R@averages,derivative(points),rtol=2e-11,atol=2e-12)
    # Arbitrary inputs retain all face traces, not just this polynomial.
    random=np.random.default_rng(4).normal(size=len(edges))
    _,traces=primitive_basis(edges,edges)
    np.testing.assert_allclose(traces@random,random,rtol=2e-12,atol=2e-12)
    np.testing.assert_allclose(R@(D@faces),derivative(points),rtol=2e-11,atol=2e-12)


def test_actual_continuity_weak_identity_and_physical_face_adjoints():
    _,model,thermal,flux=initialize();geo=model.geometry;mesh=model.mesh;rng=np.random.default_rng(25)
    radial=rng.normal(size=(mesh.nr+1,mesh.nz));radial[[0,-1]]=0
    faces=ParticleFaces(radial,rng.normal(size=(mesh.nr,mesh.nz)))
    nt=flux.density_rate(faces);particle=flux.sample(faces)
    phi=rng.normal(size=thermal.scalar.shape)
    nx,dx,_=direct_basis(.5*mesh.rf**2,.5*geo.radial_points**2,5)
    nz,dz,_=direct_basis(mesh.zf,geo.axial_points,5,True)
    samples=geo.tensor(nx,phi,nz)
    gradient=np.stack((geo.r*geo.tensor(dx,phi,nz),np.zeros_like(geo.r),geo.tensor(nx,phi,dz)),axis=-1)
    first=geo.weight*model.density_map.sample(nt)*samples
    second=geo.weight*np.sum(particle*gradient,axis=-1)
    assert abs(np.sum(first)-np.sum(second))<3e-12*np.sum(abs(first)+abs(second))
    fr,fz=flux.integrated(faces)
    independent=-(np.diff(fr,axis=0)+np.roll(fz,-1,axis=1)-fz)/mesh.volumes
    np.testing.assert_allclose(nt,independent,rtol=2e-13,atol=2e-10)
    assert abs(np.sum(mesh.volumes*nt))<1e-14*np.sum(mesh.volumes*abs(nt))
    value=rng.normal(size=particle.shape);load=flux.vector_adjoint(value)
    actual=np.sum(faces.radial*load.radial)+np.sum(faces.axial*load.axial)
    wanted=np.sum(geo.weight*np.sum(particle*value,axis=-1))
    np.testing.assert_allclose(actual,wanted,rtol=5e-13)
    q=rng.normal(size=nt.shape);adj=flux.density_adjoint(q)
    np.testing.assert_allclose(np.sum(faces.radial*adj.radial)+np.sum(faces.axial*adj.axial),np.sum(nt*q),rtol=5e-13)


def test_bulk_entropy_and_energy_work_for_both_species_with_actual_compression():
    ref,model,electron,flux=initialize();geo=model.geometry;mesh=model.mesh
    momentum=cell_average(mesh,lambda r,z:ref.fields(r,z)['particle'][...,[0,2]],order=12)
    faces=flux.from_momentum(momentum)
    for gamma in (5/3,1.4):
        thermal=CompatibleLogThermal(model,gamma=gamma,scalar=electron.scalar)
        r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');f=ref.fields(r,z)
        p=ref.thermal(r,z,gamma=gamma)['pressure']
        s=np.log(p/(f['n']*E*(f['n']/1e20)**(gamma-1)))
        rate=CompatibleBulkThermal(thermal,flux).rate(s,faces)
        energy=thermal.energy_rate(s,rate.entropy,rate.density)
        power=np.sum(faces.radial*rate.face_work.radial)+np.sum(faces.axial*rate.face_work.axial)
        scale=np.sum(mesh.volumes*(abs(energy['entropy'])+abs(energy['density'])))
        assert abs(np.sum(mesh.volumes*energy['total'])-power)<2e-12*scale
        entropy=thermal.entropy_rate(s,rate.entropy,rate.density)
        scale=np.sum(geo.weight*(abs(model.density*thermal.scalar.sample(rate.entropy))+abs(model.density_map.sample(rate.density)*thermal.scalar.sample(s))))
        assert abs(np.sum(mesh.volumes*entropy))<2e-12*scale
        constant=CompatibleBulkThermal(thermal,flux).rate(np.full(s.shape,2.),faces)
        np.testing.assert_array_equal(constant.entropy,np.zeros(s.shape))


def test_manufactured_continuity_and_bulk_thermal_rate_converge_physically():
    errors=[];axis=[]
    for nr in (8,16,32):
        ref,model,thermal,flux=initialize(nr);geo=model.geometry;mesh=model.mesh
        momentum=cell_average(mesh,lambda r,z:ref.fields(r,z)['particle'][...,[0,2]],order=12)
        faces=flux.from_momentum(momentum)
        r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');n=ref.fields(r,z)['n'];p=ref.thermal(r,z)['pressure']
        s=np.log(p/(n*E*(n/1e20)**(thermal.gamma-1)))
        rate=CompatibleBulkThermal(thermal,flux).rate(s,faces)
        f=ref.fields(geo.r,geo.z);t=ref.thermal(geo.r,geo.z)
        gradp=-t['electric']*(E*f['n'])[...,None]
        grads=gradp/t['pressure'][...,None]-thermal.gamma*f['gradn']/f['n'][...,None]
        st=-np.sum(f['particle']*grads,axis=-1)/f['n']
        exact=t['pressure']/(thermal.gamma-1)*(st+thermal.gamma*f['nt']/f['n'])
        actual=thermal.fields(s)['internal']*(thermal.scalar.sample(rate.entropy)+thermal.gamma*model.density_map.sample(rate.density)/model.density)
        def norm(a,b,part=slice(None)):
            return np.sqrt(np.sum(geo.weight[part]*(a[part]-b[part])**2)/np.sum(geo.weight[part]*b[part]**2))
        errors.append([norm(model.density_map.sample(rate.density),f['nt']),norm(thermal.scalar.sample(rate.entropy),st),norm(actual,exact)])
        axis.append(norm(thermal.scalar.sample(rate.entropy),st,slice(0,1)))
    errors=np.array(errors)
    assert np.all(errors[1]<.2*errors[0]) and np.all(errors[2]<.2*errors[1]),errors
    assert axis[1]<.2*axis[0] and axis[2]<.2*axis[1],axis


def test_boundary_quadrature_and_atomic_contracts():
    _,model,thermal,flux=initialize();shape=thermal.scalar.shape;s=np.full(shape,2.)
    before=s.copy();bad=ParticleFaces(np.ones((model.mesh.nr+1,model.mesh.nz)),np.zeros((model.mesh.nr,model.mesh.nz)))
    with pytest.raises(ValueError,match='axis/wall'):CompatibleBulkThermal(thermal,flux).rate(s,bad)
    np.testing.assert_array_equal(s,before)
    low=VariableCompatibleCurrent(model.mesh,model.density_map.cell_means,model.ion_mass)
    with pytest.raises(ValueError,match='quadrature'):CompatibleParticleFlux(low.density_map)
