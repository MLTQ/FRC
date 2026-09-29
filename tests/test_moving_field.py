"""Independent geometry, whole-space magnetic energy and virtual-force checks."""
import numpy as np
import pytest
from scipy.constants import mu_0 as MU0
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.moving_field import MovingField
from plasma_bench.moving_plasma import MovingPlasma,MovingState
from plasma_bench.waveforms import Gaussian


def test_material_sphere_volume_and_exact_curved_matching():
    errors = []
    for n in (4,8,16):
        d = MovingDomain(n,2*n)
        elements = d.elements(d.fluid_reference)
        matrix = d.stiffness(elements)
        errors.append(abs(d.volumes(elements).sum()/(4*np.pi/3)-1))
        assert np.max(abs((matrix-matrix.T).data)) < 2e-13
        if n == 4:
            assert np.linalg.eigvalsh(matrix.toarray()).min() > 0
        # Material surface is polygonal; the whole field domain is the true sphere.
        full_volume = 2*np.pi*sum(np.sum(e.weights*e.points[...,0]) for e in elements)
        assert abs(full_volume/(4*np.pi*d.matching_radius**3/3)-1) < 1e-11
    assert np.all(np.asarray(errors[:-1])/errors[1:] > 3.9)


def test_induced_uniform_sphere_energy_and_homologous_magnetic_force():
    errors = []
    for n in (4,8,16):
        d = MovingDomain(n,2*n); b = np.sqrt(MU0)
        solver = MovingPlasma(d,.5*b,pressure=0.,bias=.3*b)
        scale = .9
        state = MovingState(scale*d.fluid_reference,np.zeros_like(d.fluid_reference))
        e = solver.evaluate(state)
        # Interior uniform Bp=Bp0/scale²; exterior uniform bias plus dipole.
        exact_energy = np.pi*(1/scale-.6*scale+.09*scale**3)
        exact_derivative = np.pi*(-1/scale**2-.6+.27*scale**2)
        derivative = -np.sum(e["force"]*d.fluid_reference)
        errors.append((abs(e["magnetic"].energy/exact_energy-1),abs(derivative/exact_derivative-1)))
        assert e["magnetic"].vacuum_residual < 2e-12
        assert abs(e["magnetic"].physical_energy/e["magnetic"].energy-1) < 1e-9
    errors = np.asarray(errors)
    assert np.all(errors[:-1]/errors[1:] > 3.6)


def test_full_geometric_force_matches_independent_potential_variation():
    d = MovingDomain(5,10);b=np.sqrt(MU0)
    s = MovingPlasma(d,lambda r,z:-.75*b*(1-r*r-z*z),
                    pressure=lambda r,z:5.625*r*r*(1-r*r-z*z),bias=b)
    x = d.fluid_reference.copy();x[:,0] *= 1+.07*x[:,1];x[:,1] += .03*x[:,0]**2
    state = MovingState(x,np.zeros_like(x));e=s.evaluate(state)
    direction = np.stack((x[:,0]*(.3+x[:,1]),.2*x[:,0]**2+.1),axis=-1)
    h=1e-6
    def energy(delta):
        result = s.evaluate(MovingState(x+delta*direction,state.velocity))
        return result["internal_energy"]+result["magnetic"].energy
    measured = (energy(h)-energy(-h))/(2*h)
    assert abs(measured+np.sum(e["force"]*direction)) < 5e-8


def test_prescribed_bias_power_matches_energy_derivative():
    d = MovingDomain(6,12);b=np.sqrt(MU0)
    s=MovingPlasma(d,.2*b,pressure=.1,bias=Gaussian(b,.04,.1),max_timestep=.001)
    t,h=.025,1e-6
    e=s.evaluate(time=t)
    measured=(s.evaluate(time=t+h)["magnetic"].energy-s.evaluate(time=t-h)["magnetic"].energy)/(2*h)
    assert abs(measured-e["source_power"]) < 2e-8*max(abs(measured),1.)


def test_regular_axis_and_solenoidal_flux_on_a_deformed_cell():
    d=MovingDomain(6,12);x=d.fluid_reference.copy()
    x[:,0] *= 1+.04*x[:,1];x[:,1] += .02*x[:,0]**2
    field=MovingField(d);a=.001*(1-.2*x[:,0]**2+.1*x[:,1])
    solution=field.solve(x,a)
    # Independent line quadrature of r B dot n around one actual material cell.
    cell=d.cells[3*d.nt+2];points=x[cell];values=a[cell]
    nodes,weights=np.polynomial.legendre.leggauss(12);t=(nodes+1)/2;weights=weights/2
    result=0.
    for xi,tau,dxi,dtau in ((t,np.zeros_like(t),1.,0.),(np.ones_like(t),t,0.,1.),
                             (1-t,np.ones_like(t),-1.,0.),(np.zeros_like(t),1-t,0.,-1.)):
        shape=np.stack(((1-xi)*(1-tau),xi*(1-tau),xi*tau,(1-xi)*tau),axis=-1)
        derivative=np.stack((np.stack((-(1-tau),1-tau,tau,-tau),axis=-1),
                              np.stack((-(1-xi),-xi,xi,1-xi),axis=-1)),axis=-1)
        location=shape@points;jac=np.einsum("qnb,na->qab",derivative,points)
        gradient=np.einsum("qnb,qba,n->qa",derivative,np.linalg.inv(jac),values)
        value=shape@values;r=location[:,0]
        br,bz=-r*gradient[:,1],2*value+r*gradient[:,0]
        tangent=jac@np.array([dxi,dtau])
        result += np.sum(weights*r*(br*tangent[:,1]-bz*tangent[:,0]))
    assert abs(result) < 1e-16
    assert solution.vacuum_residual < 1e-12


def test_nonphysical_geometry_is_rejected():
    d=MovingDomain(4,8)
    x=d.fluid_reference.copy();x[d.boundary] *= 2.1
    with pytest.raises(ValueError,match="matching sphere"):
        d.elements(x)
    x=d.fluid_reference.copy();x[d.ids[2,3]]=x[d.ids[3,3]]*1.1
    with pytest.raises(ValueError,match="Jacobian"):
        d.elements(x)
