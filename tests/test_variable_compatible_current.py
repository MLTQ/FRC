"""Independent density moments, physical energy, forcing and original block solve."""
import numpy as np
import pytest
from scipy.constants import atomic_mass,elementary_charge as E,mu_0 as MU0
from scipy.sparse import bmat,csc_matrix,eye
from scipy.sparse.linalg import spsolve
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.compatible_density import CompatibleDensity
from plasma_bench.variable_compatible_current import VariableCompatibleCurrent,weak_vector
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from plasma_bench.kernel_induction import KernelInduction
from plasma_bench.compatible_charge import CompatibleCharge
from plasma_bench.compatible_pressure import CompatibleSpeciesPressure
from examples.finite_mass_manufactured import cell_average
from examples.variable_current_reference import VariableCurrentReference


def setup():
    mesh=MHDMesh(.035,.12,6,8);geo=CompatibleGeometry(mesh,quadrature_order=8)
    n=cell_average(mesh,lambda r,z:1e18*(1+.7*(r/mesh.radius)**2))
    m=VariableCompatibleCurrent(mesh,n,39.948*atomic_mass,geometry=geo)
    x=.5*mesh.corners[0]**2;X=.5*mesh.radius**2
    psi=2e-4*x*(1-x/X)**3;psi[[0,-1]]=0.
    g=cell_average(mesh,lambda r,z:.02*(1-(r/mesh.radius)**2)**2)
    return m,HarmonicMagnetic(mesh,psi,g,.005)


def test_continuum_reference_curl_faraday_and_particle_continuity():
    ref=VariableCurrentReference();r,z=np.meshgrid(np.linspace(.001,.034,7),np.linspace(-.11,.11,9),indexing='ij')
    f=ref.fields(r,z);step=1e-28
    dr={k:np.imag(v)/step for k,v in ref.fields(r+1j*step,z).items()}
    dz={k:np.imag(v)/step for k,v in ref.fields(r,z+1j*step).items()}
    curlb=np.stack((-dz['magnetic'][...,1],dz['magnetic'][...,0]-dr['magnetic'][...,2],
                    dr['magnetic'][...,1]+f['magnetic'][...,1]/r),axis=-1)
    minuscurl=np.stack((dz['electric'][...,1],dr['electric'][...,2]-dz['electric'][...,0],
                        -dr['electric'][...,1]-f['electric'][...,1]/r),axis=-1)
    np.testing.assert_allclose(curlb/MU0,f['current'],rtol=5e-13,atol=1e-11)
    np.testing.assert_allclose(minuscurl,f['magnetic_rate'],rtol=5e-13,atol=1e-12)
    nt=-dr['particle'][...,0]-f['particle'][...,0]/r-dz['particle'][...,2]
    np.testing.assert_allclose(nt,f['nt'],rtol=5e-13,atol=1e8)


def test_actual_density_means_conservation_adjoint_and_atomic_positivity():
    m,_=setup();geo=m.geometry;mesh=m.mesh;mapping=m.density_map
    exact=1e18*(1+.7*(geo.r/mesh.radius)**2)
    np.testing.assert_allclose(mapping.samples,exact,rtol=3e-14)
    np.testing.assert_allclose(m.average(mapping.samples),mapping.cell_means,rtol=3e-14)
    rng=np.random.default_rng(11);rate=rng.normal(size=mapping.cell_means.shape)*1e21
    scalar=rng.normal(size=geo.r.shape)
    np.testing.assert_allclose(np.sum(mapping.adjoint(scalar)*rate),
        np.sum(geo.weight*scalar*mapping.sample(rate)),rtol=3e-14)
    with pytest.raises(ValueError):mapping.cell_means[0,0]=0.
    # Positive means alone do not justify a positive high-order reconstruction.
    invalid=np.ones_like(mapping.cell_means);invalid[2,3]=1000
    before=invalid.copy()
    with pytest.raises(ValueError,match='no floor'):CompatibleDensity(geo,invalid)
    np.testing.assert_array_equal(before,invalid)
    with pytest.raises(ValueError,match='finite numerical range'):
        VariableCompatibleCurrent(mesh,np.full_like(invalid,1e-300),m.ion_mass,geometry=geo)
    overflowing=np.ones_like(invalid);overflowing[2,3]=1e308
    with pytest.raises(ValueError,match='finite numerical range'):
        mapping.sample(overflowing)


def test_variable_inertia_physical_integral_and_energy_directional_derivative():
    m,field=setup();geo=m.geometry;mesh=m.mesh;j=field.current(m)
    rng=np.random.default_rng(3);a=rng.normal(size=j.size);b=rng.normal(size=j.size)
    np.testing.assert_allclose(a@(m.inertia@b),b@(m.inertia@a),rtol=1e-14)
    assert a@(m.inertia@a)>0
    np.testing.assert_allclose(.5*j@(m.inertia@j),np.sum(mesh.volumes*m.current_energy(j)),rtol=1e-14)
    nt=cell_average(mesh,lambda r,z:1e18*(.2-.7*(r/mesh.radius)**2))
    jt=.3*j
    actual=m.current_energy_rate(j,jt,nt)
    np.testing.assert_allclose(np.sum(mesh.volumes*actual['density']),
        np.sum(m.density_energy_gradient(j)*nt),rtol=2e-14)
    errors=[]
    for step in (2e-2,1e-2,5e-3):
        plus=VariableCompatibleCurrent(mesh,m.density_map.cell_means+step*nt,m.ion_mass,geometry=geo)
        minus=VariableCompatibleCurrent(mesh,m.density_map.cell_means-step*nt,m.ion_mass,geometry=geo)
        difference=(plus.current_energy(j+step*jt)-minus.current_energy(j-step*jt))/(2*step)
        errors.append(np.linalg.norm(difference-actual['total'])/np.linalg.norm(actual['total']))
    assert errors[1]<.26*errors[0] and errors[2]<.26*errors[1]
    assert errors[-1]<1e-5
    # Local harmonic cross work is retained despite its zero global integral.
    bt=.3*field.packed_increment(geo)
    full=m.magnetic_energy_rate(field,bt)
    zero=HarmonicMagnetic(mesh,field.psi_increment,field.g,0.)
    cross=full-m.magnetic_energy_rate(zero,bt)
    assert np.linalg.norm(cross)>np.linalg.norm(full)*.5
    assert abs(np.sum(mesh.volumes*cross))<1e-13*np.sum(mesh.volumes*abs(cross))


def test_original_mixed_variable_density_system_and_density_energy_exchange():
    m,field=setup();geo=m.geometry;solver=KernelInduction(m,tolerance=1e-10,max_iterations=2000)
    active=solver.active;C=geo.faraday[active];H=m.hessian[active][:,active]
    rng=np.random.default_rng(17);other=m.mass@rng.normal(size=m.mass.shape[0])
    current=field.current(m);nt=.3*m.density_map.cell_means
    f=other+m.density_forcing(current,nt)
    inertia=csc_matrix(m.inertia@np.eye(m.mass.shape[0]))
    original=bmat([[eye(len(active)),-C,None],[None,m.mass,-inertia],[C.T@H,None,m.mass]],format='csc')
    direct=spsolve(original,np.r_[np.zeros(len(active)),f,np.zeros(m.mass.shape[0])])
    result=solver.solve(f)
    np.testing.assert_allclose(result.magnetic[active],direct[:len(active)],rtol=3e-9,atol=1e-8)
    assert result.residual<=1e-10
    change=m.current_energy_rate(current,result.current,nt)
    rate=np.sum(m.mesh.volumes*(m.magnetic_energy_rate(field,result.magnetic)+change['total']))
    expected=-current@other-np.sum(m.mesh.volumes*change['density'])
    assert abs(rate-expected)<5e-10*(abs(current)@abs(other)+abs(expected))


def test_manufactured_physical_polynomial_moving_density_induction():
    m,field=setup();geo=m.geometry;mesh=m.mesh;X=.5*mesh.radius**2;y=.5*geo.r**2/X
    j=field.current(m);n=1e18*(1+.7*y);nt=cell_average(mesh,lambda r,z:1e21*(1-.4*(r/mesh.radius)**2))
    nt_exact=1e21*(1-.4*y);tau=1e-5
    physical_j=np.stack((np.zeros_like(y),2e-4*geo.r*(6-18*y+12*y*y)/(X*MU0),
                         .04*(1-4*y+3*y*y)/MU0),axis=-1)
    electric=np.stack((np.zeros_like(y),-1e-4*geo.r*(1-y)**3/tau,
                       -.02*X*(1-y)**3/(3*tau)),axis=-1)
    coefficient=m.reduced_mass/(n*E**2)
    forcing=electric-coefficient[...,None]*physical_j/tau
    forcing+=(coefficient*nt_exact/n)[...,None]*physical_j
    solved=KernelInduction(m,tolerance=1e-10,max_iterations=2000).solve(
        weak_vector(geo,forcing)+m.density_forcing(j,nt))
    expected=field.packed_increment(geo)/tau
    np.testing.assert_allclose(solved.magnetic,expected,rtol=3e-8,atol=1e-9)
    np.testing.assert_allclose(geo.sample_electric(solved.current),physical_j/tau,rtol=4e-8,atol=.01)
    np.testing.assert_allclose(geo.sample_electric(solved.electric),electric,rtol=4e-8,atol=1e-10)


def test_charge_graph_uses_physical_hessian_and_preserves_local_divergence():
    m,field=setup();geo=m.geometry;mesh=m.mesh
    g=field.g*(1+.3*np.sin(2*np.pi*mesh.z[None,:]/(2*mesh.half_length)))
    charge=CompatibleCharge(geo);ir,iz=charge.current_flux(g)
    j=m.current(geo.pack(field.psi_increment,g));weak=m.mass@j
    er,_,ez=geo.split_electric(weak)
    np.testing.assert_allclose(ir[1:-1],er.reshape(mesh.nr,mesh.nz)/mesh.radial_measure[:,None],rtol=2e-12,atol=1e-14)
    np.testing.assert_allclose(iz,ez.reshape(mesh.nr+1,mesh.nz)/mesh.dz,rtol=2e-12,atol=1e-14)
    divergence=charge.divergence(ir,iz)*charge.volumes
    scale=abs(ir[:-1])+abs(ir[1:])+abs(iz)+abs(np.roll(iz,1,axis=1))
    assert np.max(abs(divergence))<5*np.finfo(float).eps*np.max(scale)


@pytest.mark.parametrize('gamma',(5/3,1.4))
def test_both_species_pressure_work_and_fixed_density_entropy(gamma):
    m,field=setup();geo=m.geometry;mesh=m.mesh
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');phase=np.pi*z/mesh.half_length
    n=1e18*(1+.3*(r/mesh.radius)**2+.2*np.sin(phase))
    bridge=CompatibleSpeciesPressure(geo,n,gamma=gamma)
    u=bridge.thermo.from_pressure(10*(1+.2*np.cos(phase)+.1*(r/mesh.radius)**2))
    g=field.g*(1+.3*np.sin(np.pi*mesh.z[None,:]/mesh.half_length))
    pair=bridge.pair(g,u);current=m.current(geo.pack(field.psi_increment,g))
    work=current@(m.mass@pair.electric_full);thermal=bridge.thermal_energy_rate(pair.entropy_rate)
    assert abs(work-thermal)<2e-12*max(abs(work),abs(thermal))
    vol=bridge.geometry.volumes
    log_rate=vol*n/u*pair.entropy_rate
    assert abs(log_rate.sum())<5e-14*np.sum(abs(log_rate))
    assert abs(np.sum(vol*pair.entropy_rate))<5e-14*np.sum(abs(vol*pair.entropy_rate))
