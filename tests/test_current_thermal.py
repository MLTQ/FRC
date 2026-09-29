"""Independent species momentum, actual finite-inertia waves and energy contracts."""
from dataclasses import replace
import numpy as np
import pytest
from scipy.constants import elementary_charge as E, electron_mass as ME, mu_0 as MU0, atomic_mass
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.current_thermal import CurrentThermal
from plasma_bench.current_thermal_midpoint import CurrentThermalMidpoint
from plasma_bench.hall_midpoint import HallConvergenceError
from examples.current_beltrami_reference import CurrentBeltrami
from examples.electron_reference import ElectronReference


@pytest.mark.parametrize('helicity', (-1, 1))
def test_finite_amplitude_reference_satisfies_both_parent_species_and_faraday(helicity):
    ref = CurrentBeltrami(helicity=helicity, amplitude=4e-5, number_density=1e19)
    r, z = np.meshgrid(np.linspace(.002, .033, 9), np.linspace(-.11, .11, 11), indexing='ij')
    f = ref.fields(r, z); alpha, beta, reduced = ref.masses; n = ref.number_density
    ion_lhs = ref.ion_mass*n*(alpha*f['drift_rate']+alpha**2*f['drift_convective'])
    electron_lhs = ME*n*(-beta*f['drift_rate']+beta**2*f['drift_convective'])
    ion_rhs = n*E*(f['electric']+np.cross(f['ion_velocity'], f['magnetic']))+beta*f['reaction']
    electron_rhs = -n*E*(f['electric']+np.cross(f['electron_velocity'], f['magnetic']))+alpha*f['reaction']
    for a, b in ((ion_lhs, ion_rhs), (electron_lhs, electron_rhs)):
        assert np.linalg.norm(a-b) < 3e-12*np.linalg.norm(a)
    power = np.sum(f['reaction']*(beta*f['ion_velocity']+alpha*f['electron_velocity']), axis=-1)
    assert np.linalg.norm(power) < 1e-13*np.linalg.norm(f['reaction'])*np.linalg.norm(f['drift'])
    h = 1e-5*ref.radius
    def derivative(coordinate):
        values = [ref.fields(r+i*h if coordinate == 0 else r, z+i*h if coordinate == 1 else z)['electric'] for i in (-2,-1,1,2)]
        return (values[0]-8*values[1]+8*values[2]-values[3])/(12*h)
    dr, dz = derivative(0), derivative(1)
    minus_curl = np.stack((dz[..., 1], dr[..., 2]-dz[..., 0], -dr[..., 1]-f['electric'][..., 1]/r), axis=-1)
    assert np.linalg.norm(minus_curl-f['magnetic_rate']) < 2e-8*np.linalg.norm(f['magnetic_rate'])
    kinetic = .5*n*(ref.ion_mass*np.sum(f['ion_velocity']**2,axis=-1)+ME*np.sum(f['electron_velocity']**2,axis=-1))
    np.testing.assert_allclose(kinetic, .5*n*reduced*np.sum(f['drift']**2,axis=-1), rtol=5e-16)
    def total_flux(rr, zz):
        ff = ref.fields(rr, zz)
        kk = .5*n*reduced*np.sum(ff['drift']**2, axis=-1)
        enthalpy = 2.5*n*E*(alpha*ref.ion_ev-beta*ref.electron_ev)
        return np.cross(ff['electric'], ff['magnetic'])/MU0+ff['drift']*(enthalpy-(beta-alpha)*kk)[..., None]
    def flux_derivative(axis):
        values = [total_flux(r+i*h if axis == 0 else r, z+i*h if axis == 1 else z) for i in (-2,-1,1,2)]
        return (values[0]-8*values[1]+8*values[2]-values[3])/(12*h)
    energy_rate = np.sum(f['magnetic']*f['magnetic_rate'],axis=-1)/MU0+n*reduced*np.sum(f['drift']*f['drift_rate'],axis=-1)
    divergence = flux_derivative(0)[..., 0]+total_flux(r,z)[..., 0]/r+flux_derivative(1)[..., 2]
    assert np.linalg.norm(energy_rate+divergence) < 2e-7*np.linalg.norm(energy_rate)


def initialize(ref, nr):
    mesh = MHDMesh(ref.radius, ref.length/2, nr, 2*nr)
    n = np.full((nr+1,2*nr), ref.number_density)
    model = CurrentThermal(mesh,n,ref.ion_mass,model=ref.model)
    ue = model.electron.from_pressure(n*ref.electron_ev*E); ui = model.ion.from_pressure(n*ref.ion_ev*E)
    return model, *ref.initial_fields(mesh), ue, ui


def test_actual_finite_mass_nonlinear_wave_refines_with_energy_and_two_entropies():
    ref = CurrentBeltrami(number_density=1e18, amplitude=1e-7)
    errors=[]
    for nr in (4,8,16):
        model,psi,g,ue,ui = initialize(ref,nr); geo=model.geometry
        initial=model.diagnostics(geo.pack(psi,g),ue,ui)
        solver=CurrentThermalMidpoint(model); dt=.25/abs(ref.omega)/20
        flux=np.zeros_like(g)
        for _ in range(20):
            result=solver.step(psi,g,ue,ui,dt)
            psi,g,ue,ui=result.psi,result.g,result.electron_entropy,result.ion_entropy
            flux+=result.total_energy_change
        final=model.diagnostics(geo.pack(psi,g),ue,ui)
        exact=ref.initial_fields(model.mesh,20*dt)[1]
        errors.append(np.sqrt(np.sum(model.mesh.volumes*(g-exact)**2)/np.sum(model.mesh.volumes*exact**2)))
        assert abs(final['total_energy_j']/initial['total_energy_j']-1)<2e-12
        assert abs(np.sum(model.mesh.volumes*flux))<1e-15*initial['total_energy_j']
        for species in ('electron','ion'):
            assert final['min_'+species+'_pressure_pa']>0
            assert abs(final[species+'_entropy_proxy_j']/initial[species+'_entropy_proxy_j']-1)<2e-14
    assert errors[1]<.6*errors[0] and errors[2]<.6*errors[1],errors


def test_nonlinear_variable_density_joint_work_and_atomic_rejection():
    ref=ElectronReference();mesh=MHDMesh(ref.radius,ref.length/2,8,16)
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');f=ref.fields(r,z)
    model=CurrentThermal(mesh,f['density'],39.948*atomic_mass)
    psi=np.zeros((9,17));g=ref.cell_average(mesh,'g')
    ue=model.electron.from_pressure(f['pressure']);ui=model.ion.from_pressure(np.full_like(ue,16.))
    arrays=[v.copy() for v in (psi,g,ue,ui)]
    before=model.diagnostics(model.geometry.pack(psi,g),ue,ui)
    result=CurrentThermalMidpoint(model).step(psi,g,ue,ui,1e-8)
    after=model.diagnostics(model.geometry.pack(result.psi,result.g),result.electron_entropy,result.ion_entropy)
    assert abs(after['total_energy_j']/before['total_energy_j']-1)<2e-13
    assert abs(after['current_kinetic_energy_j']-before['current_kinetic_energy_j'])>0
    with pytest.raises(HallConvergenceError):
        CurrentThermalMidpoint(model,tolerance=1e-13,max_iterations=1).step(psi,g,ue,ui,1e-7)
    for old,actual in zip(arrays,(psi,g,ue,ui)):
        np.testing.assert_array_equal(old,actual)


def test_exact_massless_and_pair_mass_coefficients_are_not_effective_mass_fits():
    mesh=MHDMesh(.035,.12,4,8);n=np.full((5,8),1e20)
    hall=CurrentThermal(mesh,n,39.948*atomic_mass,model='hall')
    pair=CurrentThermal(mesh,n,ME)
    assert hall.alpha==0 and hall.beta==1 and hall.delta==1 and hall.reduced_mass==0 and hall.inertia.nnz==0
    assert pair.alpha==.5 and pair.beta==.5 and pair.delta==0 and pair.reduced_mass==ME/2
    with pytest.raises(ValueError):
        hall.density[0,0,0,0]=0


def test_independent_physical_current_and_species_energy_quadrature():
    from examples.current_thermal_norms import physical_energy
    mesh=MHDMesh(.035,.12,8,16)
    f=ElectronReference().fields(*np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij'))
    model=CurrentThermal(mesh,f['density'],39.948*atomic_mass)
    psi,g=CurrentBeltrami(amplitude=4e-5).initial_fields(mesh)
    ue=model.electron.from_pressure(f['pressure']);ui=model.ion.from_pressure(.2*f['pressure'])
    packed=model.geometry.pack(psi,g);actual=model.energy_density(packed,ue,ui)
    independent=[]
    for order in (8,14):
        values=physical_energy(mesh,psi,g,model.current(packed),model.electron.density,ue,ui,
            reduced_mass=model.reduced_mass,order=order)
        independent.append(values)
        for key in values:
            assert np.linalg.norm(values[key]-actual[key])<2e-10*np.linalg.norm(values[key])
    for key in actual:
        assert np.linalg.norm(independent[0][key]-independent[1][key])<2e-13*np.linalg.norm(actual[key])
