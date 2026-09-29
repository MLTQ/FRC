"""Independent polynomial, integral, physical-current and mixed-equation checks."""
import numpy as np
from scipy.constants import atomic_mass,mu_0 as MU0
from scipy.sparse import bmat,eye
from scipy.sparse.linalg import spsolve
from plasma_bench.compatible_polynomials import bases
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.compatible_current import CompatibleCurrent
from plasma_bench.compatible_velocity import velocity
from plasma_bench.kernel_induction import KernelInduction
from plasma_bench.mhd_mesh import MHDMesh
from examples.finite_mass_manufactured import cell_average


def test_interval_integrals_and_independent_polynomial_derivatives():
    for degree in (1,3,5):
        edges=np.linspace(0.,1.,9)**2;q,w=np.polynomial.legendre.leggauss(9)
        points=(.5*(edges[:-1]+edges[1:])[:,None]+.5*np.diff(edges)[:,None]*q).ravel()
        N,dN,E,D=bases(edges,points,degree)
        for power in range(degree+1):
            values=edges**power
            np.testing.assert_allclose(N@values,points**power,rtol=2e-12,atol=2e-13)
            target=power*points**(power-1) if power else np.zeros_like(points)
            np.testing.assert_allclose(dN@values,target,rtol=3e-11,atol=2e-11)
        integrated=np.einsum('iqk,q->ik',E.toarray().reshape(8,9,8),w)/2
        np.testing.assert_allclose(integrated,np.eye(8),rtol=0.,atol=2e-12)
        np.testing.assert_allclose((E@D).toarray(),dN.toarray(),rtol=2e-12,atol=2e-11)
    edges=np.linspace(-1.,1.,13)
    points=(.5*(edges[:-1]+edges[1:])[:,None]+.5*np.diff(edges)[:,None]*q).ravel()
    N,dN,E,D=bases(edges,points,5,periodic=True)
    integrated=np.einsum('iqk,q->ik',E.toarray().reshape(12,9,12),w)/2
    np.testing.assert_allclose(integrated,np.eye(12),rtol=0.,atol=3e-13)


def test_independent_polynomial_magnetic_field_current_and_energy():
    mesh=MHDMesh(.7,.6,8,12);m=CompatibleCurrent(mesh,1e18,39.948*atomic_mass);geo=m.geometry
    x=.5*mesh.corners[0]**2;X=.5*mesh.radius**2;y=x/X
    psi=.4*x+.7*x*(1-y)**3
    g=cell_average(mesh,lambda r,z:.3*(1-.5*r*r/X)**2)
    magnetic=geo.pack(psi,g);current=m.current(magnetic)
    y=.5*geo.r**2/X
    exact_b=np.stack((np.zeros_like(y),.3*geo.r*(1-y)**2,.4+.7*(1-6*y+9*y*y-4*y**3)),axis=-1)
    exact_j=np.stack((np.zeros_like(y),.7*geo.r*(6-18*y+12*y*y)/(X*MU0),.6*(1-4*y+3*y*y)/MU0),axis=-1)
    np.testing.assert_allclose(geo.magnetic_samples(psi,g),exact_b,rtol=3e-11,atol=1e-12)
    np.testing.assert_allclose(geo.sample_electric(current),exact_j,rtol=2e-9,atol=2e-6)
    exact_energy=np.sum(geo.weight*np.sum(exact_b**2,axis=-1))/(2*MU0)
    np.testing.assert_allclose(.5*magnetic@(m.hessian@magnetic),exact_energy,rtol=1e-13)
    rng=np.random.default_rng(42);electric=rng.normal(size=m.mass.shape[0]);rate=geo.faraday@electric
    work=magnetic@(m.hessian@rate)+current@(m.mass@electric)
    scale=abs(magnetic)@abs(m.hessian@rate)+abs(current)@abs(m.mass@electric)
    assert abs(work)<50*np.finfo(float).eps*scale


def test_original_faraday_polynomial_and_curl_free_potential():
    mesh=MHDMesh(.7,.6,8,12);geo=CompatibleGeometry(mesh);X=.5*mesh.radius**2
    x=.5*mesh.rf**2
    # Ez=1+2x+3x², Etheta=r(1-x/X), Er=r; all are exact polynomial fields.
    er=np.ones((mesh.nr,mesh.nz));et=np.broadcast_to((mesh.rf[1:-1]*(1-x[1:-1]/X))[:,None],(mesh.nr-1,mesh.nz))
    ez=np.broadcast_to((1+2*x+3*x*x)[:,None],(mesh.nr+1,mesh.nz))
    electric=np.r_[er.ravel(),et.ravel(),ez.ravel()]
    rate=geo.magnetic_samples(*geo.unpack(geo.faraday@electric));x=.5*geo.r**2
    exact=np.stack((np.zeros_like(x),geo.r*(2+6*x),-2+4*x/X),axis=-1)
    np.testing.assert_allclose(rate,exact,rtol=2e-12,atol=5e-12)
    phi=(1+(.5*mesh.rf[:,None]**2)**3)*np.cos(2*np.pi*mesh.zf[None,:-1]/(2*mesh.half_length))
    gradient=geo.gradient(phi);defect=geo.faraday@gradient
    scale=abs(geo.faraday)@abs(gradient)
    assert np.linalg.norm(defect)<2e-15*np.linalg.norm(scale)


def test_conservative_variable_density_velocity_polynomial():
    mesh=MHDMesh(.7,.6,8,12);geo=CompatibleGeometry(mesh)
    def exact(r,z):
        rho=2+.5*r*r
        return np.stack((rho,rho*r*(1+r*r),rho*(3+r*r),np.zeros_like(r),rho*r*r*(2+r*r)),axis=-1)
    state=cell_average(mesh,exact)
    target=np.stack((geo.r*(1+geo.r**2),geo.r*(2+geo.r**2),3+geo.r**2),axis=-1)
    np.testing.assert_allclose(velocity(geo,state),target,rtol=3e-13,atol=2e-13)


def test_independent_original_three_block_induction():
    mesh=MHDMesh(.035,.12,6,8);m=CompatibleCurrent(mesh,1e18,39.948*atomic_mass,quadrature_order=8)
    solver=KernelInduction(m,tolerance=1e-10,max_iterations=2000);active=solver.active
    C=m.geometry.faraday[active];H=m.hessian[active][:,active]
    rng=np.random.default_rng(4);forcing=m.mass@rng.normal(size=m.mass.shape[0])
    original=bmat([[eye(len(active)), -C,None],[None,m.mass,-m.inertia],[C.T@H,None,m.mass]],format='csc')
    direct=spsolve(original,np.r_[np.zeros(len(active)),forcing,np.zeros(m.mass.shape[0])])
    solved=solver.solve(forcing)
    np.testing.assert_allclose(solved.magnetic[active],direct[:len(active)],rtol=2e-9,atol=1e-10)
    assert solved.residual<=1e-10
