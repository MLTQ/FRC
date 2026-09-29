"""Independent magnetic quadrature, natural trace and work checks."""
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.geometry import Grid
from plasma_bench.mhd_state import from_primitive, primitive, NonphysicalState
from plasma_bench.magnetic_quadrature import energy_density, stiffness_action
from plasma_bench.coupled_benchmarks import aligned_meshes
from plasma_bench.coupling_mesh import CouplingMesh
from plasma_bench.vacuum_energy import variational_trace
from plasma_bench.vacuum_flux import flux_components
from plasma_bench.conservative_benchmarks import from_split_initial
from plasma_bench.shock_benchmarks import standing_shock
from plasma_bench.flux_limiting import admissible_update


def test_uniform_field_energy_and_independent_polynomial_integral():
    errors = []
    for n in (8,16,32):
        mesh = MHDMesh(1.2,.7,n,2*n)
        r,z = mesh.corners
        field = np.sqrt(MU0)
        uniform = .5*field*r*r+0*z
        exact_uniform = np.pi/MU0*field**2*mesh.radius**2*mesh.half_length
        assert_allclose(np.sum(energy_density(mesh,uniform)*mesh.volumes), exact_uniform, rtol=2e-15)
        # psi=a r² z gives Br=-a r, Bz=2 a z. Integrate B² r analytically.
        a = .2*np.sqrt(MU0)
        exact = np.pi*a*a/MU0*(.5*mesh.radius**4*mesh.half_length+
                               4/3*mesh.radius**2*mesh.half_length**3)
        actual = np.sum(energy_density(mesh,a*r*r*z)*mesh.volumes)
        errors.append(abs(actual-exact))
    assert_allclose(np.array(errors[:-1])/errors[1:], 4., rtol=2e-12)


def test_positive_quadrature_matches_vacuum_stiffness_for_compact_field():
    mesh,grid = aligned_meshes(12)
    mapping = CouplingMesh(mesh,grid)
    r,z = mesh.corners
    psi = np.sqrt(MU0)*r*r*(1-r*r)**3*(1-z*z)**3
    full = np.zeros(grid.shape)
    full[:mesh.nr+1,mapping.z_offset:mapping.z_offset+mesh.nz+1] = psi
    u = mapping.interior(full)
    assert energy_density(mesh,psi).min() >= 0
    assert_allclose(np.sum(energy_density(mesh,psi)*mesh.volumes), np.pi/MU0*(u@(mapping.k@u)), rtol=2e-14)
    direction = np.sqrt(MU0)*r*r*(.3+np.cos(2*z))
    epsilon = 1e-3
    fd = np.sum((energy_density(mesh,psi+epsilon*direction)-
                 energy_density(mesh,psi-epsilon*direction))*mesh.volumes)/(2*epsilon)
    analytic = 2*np.pi/MU0*np.sum(stiffness_action(mesh,psi)*direction)
    assert_allclose(fd,analytic,rtol=2e-12)


@pytest.mark.parametrize("axial_spacing_ratio", [1., .5])
def test_variational_traces_converge_to_independent_harmonic_fields_including_corners(axial_spacing_ratio):
    errors = []
    a,b,c = .7,.3,.05
    for n in (8,16,32):
        mesh,grid = aligned_meshes(n)
        if axial_spacing_ratio == .5:
            mesh,grid = MHDMesh(1.,1.,n,4*n),Grid(2.,2.,2*n+1,8*n+1)
        mapping = CouplingMesh(mesh,grid)
        r,z = grid.r[:,None],grid.z[None,:]
        psi = a*r*r+b*r*r*z+c*r*r*(r*r-4*z*z)
        load = mapping.k@mapping.grid_vector(psi)
        trace = variational_trace(mapping,psi,load)
        exact_wall = 2*a+2*b*mesh.zf+c*(4*mesh.radius**2-8*mesh.zf**2)
        exact_low = mesh.rf*(-b-8*c*mesh.half_length)
        exact_high = mesh.rf*(-b+8*c*mesh.half_length)
        errors.append(max(abs(trace.bz_wall-exact_wall).max(),
                          abs(trace.br_bottom-exact_low).max(),abs(trace.br_top-exact_high).max()))
        assert trace.br_bottom[0] == trace.br_top[0] == 0
    assert errors[1] < errors[0]/3.5 and errors[2] < errors[1]/3.5


def test_uniform_field_prescribed_electric_curl_changes_magnetic_not_thermal_energy():
    mesh,grid = aligned_meshes(8)
    mapping = CouplingMesh(mesh,grid)
    field,alpha = np.sqrt(MU0),.01*np.sqrt(MU0)
    full = .5*field*grid.r[:,None]**2+np.zeros((1,grid.nz))
    psi = mapping.fluid_flux(full)
    q = np.zeros((mesh.nr,mesh.nz,6));q[...,0]=q[...,3]=1.;q[...,5]=field
    electric = alpha*mesh.rf[:,None]*np.cos(np.pi*mesh.zf[None,:]/2)
    trace = variational_trace(mapping,full,np.zeros(mapping.node_count))
    def rhs(e):
        flux = flux_components(mesh,q,psi,5/3,e,trace.face_fields(mesh),np.zeros_like(psi))
        flux.boundary_poynting = trace
        return flux.rhs(mesh)
    du,dpsi,pout = rhs(electric)
    base,bpsi,bout = rhs(np.zeros_like(psi))
    du,dpsi,pout = du-base,dpsi-bpsi,pout-bout
    # At uniform Bz, the staggered magnetic-energy derivative uses the mean
    # of the two axial-face Bz rates. J=0, so instantaneous thermal rate is zero.
    dbz = mesh.face_fields(dpsi)[1]
    mdot = field*.5*(dbz[:,:-1]+dbz[:,1:])/MU0
    assert_allclose(du[...,3]-mdot,0,atol=3e-15)
    assert_allclose(np.sum(mdot*mesh.volumes),-pout,atol=3e-15)


def test_energy_conventions_are_explicit_and_do_not_repair_evolved_states():
    mesh = MHDMesh(1.,1.,16,16)
    r,z = mesh.corners
    i,j = np.indices(r.shape)
    psi = .03*np.sqrt(MU0)*r*r*(-1.)**(i+j)
    legacy = from_primitive(mesh,psi,1.,.01)
    before = legacy.conserved.copy()
    with pytest.raises(NonphysicalState):
        primitive(mesh,legacy,5/3,energy_quadrature="staggered")
    assert_allclose(legacy.conserved,before,rtol=0,atol=0)
    constructed = from_primitive(mesh,psi,1.,.01,energy_quadrature="staggered")
    q = primitive(mesh,constructed,5/3,energy_quadrature="staggered")
    assert_allclose(q[...,3],.01,atol=2e-16)


def test_actual_shock_fluxes_satisfy_semidiscrete_response_work_without_repair():
    solver = from_split_initial(standing_shock(24)[0], energy_quadrature="staggered",
                                boundary_poynting="variational", time_order=3)
    solver.advance(.01)
    m, mesh = solver.mapping, solver.mesh
    total = m.fluid_flux(solver.total_flux())
    ku = m.k@m.interior(solver.state.induced_psi)
    before = solver.state.copy()
    for order in (1, 2):
        flux, powers = solver._fluxes(solver.state, solver.time, order)
        du, dpsi, pout = flux.rhs(mesh)
        magnetic_rate = 2*np.pi/MU0*np.sum(stiffness_action(mesh,total)*dpsi)
        nodal_rate = 2*np.pi/MU0*np.sum(m.fluid_flux(m.full(ku))*dpsi)
        assert abs(-pout-magnetic_rate+nodal_rate) < 2e-11
        # Independently difference actual response energies after exact vacuum
        # recompletion. Constant bias makes the source power exactly zero.
        assert powers[0] == 0
        epsilon = 1e-6
        response = []
        for sign in (-1, 1):
            state = solver._new_state(before.conserved+sign*epsilon*du,
                                      total+sign*epsilon*dpsi, solver.time)
            response.append(solver._integrals(state,solver.time)["response_energy_j"])
        assert abs((response[1]-response[0])/(2*epsilon)) < 5e-8
    assert_allclose(solver.state.conserved, before.conserved, rtol=0, atol=0)
    assert_allclose(solver.state.induced_psi, before.induced_psi, rtol=0, atol=0)


def test_boundary_fallback_recomputes_work_from_selected_corners():
    mesh,grid = aligned_meshes(8)
    mapping = CouplingMesh(mesh,grid)
    field = np.sqrt(MU0)
    full = .5*field*grid.r[:,None]**2+np.zeros((1,grid.nz))
    psi = mapping.fluid_flux(full)
    state = from_primitive(mesh,psi,1.,1.,energy_quadrature="staggered")
    q = primitive(mesh,state,5/3,energy_quadrature="staggered")
    trace = variational_trace(mapping,full,np.zeros(mapping.node_count))
    electric = .001*field*mesh.rf[:,None]*np.ones((1,mesh.nz+1))
    high = flux_components(mesh,q,psi,5/3,electric,trace.face_fields(mesh),np.zeros_like(psi))
    low = flux_components(mesh,q,psi,5/3,np.zeros_like(psi),trace.face_fields(mesh),np.zeros_like(psi),1)
    high.boundary_poynting = low.boundary_poynting = trace
    high.radial_source[-1,mesh.nz//2] += 1000.
    dt = .01
    rhs,report = admissible_update(mesh,state.conserved,psi,dt,high,lambda:low,energy_quadrature="staggered")
    mixed,_ = high.local_fallback(low,report.troubled)
    assert report.troubled[-1,mesh.nz//2]
    assert report.affected.mean() < .2
    assert_allclose(rhs[0][report.troubled],low.rhs(mesh)[0][report.troubled],rtol=0,atol=0)
    # Independent nodal dual integration, using the actual mixed E. Uniform
    # B has zero end work; endpoint wall nodes carry half axial measure.
    weights = np.full(mesh.nz+1,mesh.dz);weights[[0,-1]] *= .5
    exact_outward = 2*np.pi*mesh.radius*field/MU0*np.dot(weights,mixed.electric[-1])
    assert_allclose(rhs[2],exact_outward,rtol=2e-14)
    assert_allclose(np.sum(rhs[0][...,3]*mesh.volumes),-exact_outward,atol=2e-14)
