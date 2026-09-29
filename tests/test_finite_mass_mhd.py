"""Coupled moving manufactured rates, exact uniform flow and conservative steps."""
import numpy as np
import pytest
from scipy.constants import elementary_charge as E,atomic_mass
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.finite_mass_state import FiniteMassParameters,from_primitive
from plasma_bench.finite_mass_rhs import rhs,pressure_rates,FiniteMassSources
from plasma_bench.finite_mass_mhd import FiniteMassMHD
from plasma_bench.mhd_state import NonphysicalState
from examples.moving_species_reference import MovingSpeciesReference
from examples.finite_mass_manufactured import initialize,ManufacturedSources,point_conserved,cell_average


def test_coupled_manufactured_mass_momentum_energy_and_species_rates_refine():
    ref=MovingSpeciesReference();t=.23/ref.frequency;errors=[]
    for nr in (8,16,32):
        mesh=MHDMesh(ref.radius,ref.length/2,nr,2*nr);state,par=initialize(mesh,ref,t)
        rate=rhs(mesh,state,par,time=t,source=ManufacturedSources(mesh,ref))
        exact=cell_average(mesh,lambda r,z:point_conserved(ref,r,z,t,rate=True))
        relative=np.sqrt(np.sum(mesh.volumes[...,None]*(rate.state.conserved-exact)**2,axis=(0,1))/np.sum(mesh.volumes[...,None]*exact**2,axis=(0,1)))
        pedot,pidot=pressure_rates(state,rate,par)
        piref=cell_average(mesh,lambda r,z:ref.derivatives(r,z,t)[2]['pi'])
        peref=ref.derivatives(*np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij'),t)[2]['pe']
        errors.append(np.r_[relative,np.linalg.norm(pedot-peref)/np.linalg.norm(peref),
            np.sqrt(np.sum(mesh.volumes*(pidot-piref)**2)/np.sum(mesh.volumes*piref**2))])
        # Same actual mass flux, including numerical diffusion, supplies n_t.
        mapping=rate.evaluation.model.mapping
        particle=mapping.fluxes(rate.radial_flux[...,0]/par.masses.total,rate.axial_flux[...,0]/par.masses.total)
        mismatch=rate.number_density_rate+mapping.geometry.divergence(*particle)
        assert np.linalg.norm(mismatch)<2e-14*np.linalg.norm(rate.number_density_rate)
    errors=np.array(errors)
    assert np.all(errors[1:]<.65*errors[:-1]),errors


def test_uniform_translation_and_bias_are_preserved_without_thermal_repair():
    mesh=MHDMesh(.035,.12,8,16);par=FiniteMassParameters(39.948*atomic_mass)
    r,_=mesh.corners;psi=.005*r*r/2
    state=from_primitive(mesh,par,psi,0.,par.masses.total*1e20,np.full((9,16),5e20*E),1e20*E,vz=1e4)
    model=FiniteMassMHD(mesh,state,par,max_timestep=1e-10);model.advance(1e-9)
    np.testing.assert_allclose(model.state.conserved,state.conserved,rtol=2e-14,atol=1e-18)
    np.testing.assert_allclose(model.state.psi,state.psi,rtol=2e-13,atol=1e-20)
    assert np.max(abs(model.state.g))<1e-14
    np.testing.assert_allclose(model.state.electron_entropy,state.electron_entropy,rtol=2e-14)


def test_unforced_budget_and_failed_coupled_step_are_atomic():
    ref=MovingSpeciesReference();mesh=MHDMesh(ref.radius,ref.length/2,8,16)
    state,par=initialize(mesh,ref,.23/ref.frequency);model=FiniteMassMHD(mesh,state,par,max_timestep=1e-10)
    model.advance(5e-10);integrals=model.integrals()
    for key,value in model.initial_integrals.items():assert abs(integrals[key]-value)<2e-14*abs(value)
    saved=model.state.copy();diagnostics=model.diagnostics()
    def excessive_sink(time,ev):
        fluid=np.zeros_like(state.conserved);fluid[...,3]=-1e14
        return FiniteMassSources(fluid,np.zeros_like(state.electron_entropy),np.zeros(ev.model.geometry.r.shape+(3,)))
    model.source=excessive_sink
    with pytest.raises(NonphysicalState):model.step()
    for a,b in zip(saved.arrays(),model.state.arrays()):np.testing.assert_array_equal(a,b)
    assert model.diagnostics()==diagnostics
