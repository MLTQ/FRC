"""Independent moment, contact and compressible-wave checks of smooth transport."""
import numpy as np
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.dual_mass import DualMass
from plasma_bench.annular_cell_polynomial import radial_cell_matrix,axial_cell_samples
from plasma_bench.smooth_bulk_entropy import entropy_rate
from plasma_bench.smooth_finite_mass_rhs import rhs
from plasma_bench.finite_mass_rhs import pressure_rates
from examples.finite_mass_rate_measurements import protocol,case_state
from examples.finite_mass_rate_norms import measure
from examples.finite_mass_contact_reference import FiniteMassContact
from examples.finite_mass_contact_case import initialize,averages


def test_contact_satisfies_unforced_parent_species_and_total_energy():
    from scipy.constants import elementary_charge as E
    ref=FiniteMassContact();z=np.linspace(-.1,.1,29);t=.13*ref.period;h=1e-4*ref.length
    f=ref.fields(z,t);dt=ref.fields(z,t,rate=True)
    a,b,c,d=[ref.fields(z+j*h,t) for j in (-2,-1,1,2)]
    derivative={key:(a[key]-8*b[key]+8*c[key]-d[key])/(12*h) for key in f}
    for key in ('density','electron_pressure','ion_pressure','momentum','energy'):
        scale=np.linalg.norm(dt[key])+np.linalg.norm(ref.velocity*derivative[key])
        assert np.linalg.norm(dt[key]+ref.velocity*derivative[key])/scale<2e-11
    for name,charge in (('ion',E),('electron',-E)):
        force=charge*f['number_density']*f['electric_z'];gradient=derivative[name+'_pressure']
        assert np.linalg.norm(force-gradient)/np.linalg.norm(gradient)<2e-11


def test_smooth_uniform_translation_and_rejected_trial_are_atomic():
    import pytest
    from plasma_bench.finite_mass_state import FiniteMassParameters,from_primitive
    from plasma_bench.smooth_finite_mass_mhd import SmoothFiniteMassMHD
    from plasma_bench.finite_mass_rhs import FiniteMassSources
    from plasma_bench.mhd_state import NonphysicalState
    ref=FiniteMassContact();mesh=MHDMesh(ref.radius,ref.length/2,8,16);par=FiniteMassParameters(ref.ion_mass)
    state=from_primitive(mesh,par,.5*ref.bias*mesh.corners[0]**2,0.,par.masses.total*ref.number_density,
        np.full((9,16),80.),16.,vz=ref.velocity)
    model=SmoothFiniteMassMHD(mesh,state,par,max_timestep=1e-9);model.advance(3e-9)
    np.testing.assert_allclose(model.state.conserved,state.conserved,rtol=3e-14,atol=1e-18)
    for key,value in model.initial_integrals.items():assert abs(model.integrals()[key]-value)<1e-14*max(abs(value),1e-20)
    saved=model.state.copy();before=model.diagnostics()
    def sink(time,ev):
        fluid=np.zeros_like(state.conserved);fluid[...,3]=-1e15
        return FiniteMassSources(fluid,np.zeros_like(state.electron_entropy),np.zeros(ev.model.geometry.r.shape+(3,)))
    model.source=sink
    with pytest.raises(NonphysicalState):model.step()
    for x,y in zip(model.state.arrays(),saved.arrays()):np.testing.assert_array_equal(x,y)
    assert model.diagnostics()==before


def test_transverse_polynomials_preserve_arbitrary_actual_cell_means():
    mesh=MHDMesh(.035,.12,16,24);random=np.random.default_rng(581).normal(size=(mesh.nr,mesh.nz))
    nodes,weights=np.polynomial.legendre.leggauss(5)
    for power in (0,1,2):
        total=np.zeros_like(random)
        for x,w in zip(nodes,weights):
            r=mesh.r[:,None]+.5*mesh.dr*x
            total+=.5*w*r*mesh.dr/mesh.radial_measure[:,None]*(radial_cell_matrix(mesh,.5*(x+1),power=power)@random)
        np.testing.assert_allclose(total,random,rtol=1e-12,atol=3e-13)
    total=sum(.5*w*axial_cell_samples(random,.5*(x+1)) for x,w in zip(nodes,weights))
    np.testing.assert_allclose(total,random,rtol=2e-14,atol=3e-15)


def test_smooth_entropy_matches_independent_polynomial_flux_on_dual_volumes():
    mesh=MHDMesh(.035,.12,16,24);mapping=DualMass(mesh);edges=mapping.geometry.dual_x_faces
    lo,hi=edges[:-1]/edges[-1],edges[1:]/edges[-1]
    mean=lambda degree:(hi**(degree+1)-lo**(degree+1))/((degree+1)*(hi-lo))
    n=np.broadcast_to((2+.2*mean(1)+.1*mean(2))[:,None],mapping.volumes.shape)
    u=np.broadcast_to((3-.3*mean(1)+.4*mean(3))[:,None],mapping.volumes.shape)
    x=edges/edges[-1];flow=np.broadcast_to((x*(1-x))[:,None],(mesh.nr+2,mesh.nz)).copy()
    axial=np.zeros_like(n)
    expected=-(flow[1:]*(3-.3*x[1:,None]+.4*x[1:,None]**3)/(2+.2*x[1:,None]+.1*x[1:,None]**2)
               -flow[:-1]*(3-.3*x[:-1,None]+.4*x[:-1,None]**3)/(2+.2*x[:-1,None]+.1*x[:-1,None]**2))/mapping.volumes
    np.testing.assert_allclose(entropy_rate(mapping,n,u,flow,axial),expected,rtol=2e-12,atol=1e-7)
    np.testing.assert_allclose(entropy_rate(mapping,n,7*n,flow,axial),-7*mapping.geometry.divergence(flow,axial),rtol=2e-12,atol=1e-7)


def test_shared_transport_contact_pressure_rates_refine_without_background_masking():
    ref=FiniteMassContact();errors=[]
    for nz in (16,32,64):
        mesh=MHDMesh(ref.radius,ref.length/2,8,nz);state,par=initialize(mesh,ref)
        result=rhs(mesh,state,par);pe,pi=pressure_rates(state,result,par)
        exact=averages(mesh,ref,rate=True);pereference=np.broadcast_to(ref.fields(mesh.zf[:-1],rate=True)['electron_pressure'][None],pe.shape)
        errors.append([np.linalg.norm(pe-pereference)/np.linalg.norm(pereference),np.linalg.norm(pi-exact['ion_pressure'])/np.linalg.norm(exact['ion_pressure'])])
    errors=np.array(errors)
    assert np.all(errors[1:]<.35*errors[:-1]),errors


def test_slow_and_fast_pressure_rates_refine_against_parent_species():
    lookup={c['name']:c for c in protocol()}
    for branch in (0,2):
        errors=[]
        for nr in (8,16,32):
            ref,mesh,state,par,source,t=case_state(lookup[f'wave_{branch}_{nr}'])
            result=rhs(mesh,state,par);metrics,_=measure(mesh,ref,state,par,result,t,family='wave')
            errors.append([metrics[k]['relative_L2'] for k in ('electron_pressure','ion_pressure')])
        errors=np.array(errors)
        assert np.all(errors[1:]<.4*errors[:-1]),errors
        assert max(errors[-1])<.1,errors
