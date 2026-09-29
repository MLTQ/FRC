"""Exact electric kernels and weak-signal independent mixed-system comparison."""
import numpy as np
from scipy.constants import atomic_mass
from scipy.sparse import bmat
from scipy.sparse.linalg import splu
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.current_thermal import CurrentThermal
from plasma_bench.inertial_induction import InertialInduction
from plasma_bench.kernel_induction import KernelInduction,difference_curl


def test_difference_curl_annihilates_analytic_kernel_and_matches_matrix():
    mesh=MHDMesh(.035,.12,6,12);model=CurrentThermal(mesh,np.full((7,12),1e18),39.948*atomic_mass);g=model.geometry
    radial=np.broadcast_to((1+mesh.r**2)[:,None],(mesh.nr,mesh.nz));axial=np.broadcast_to(np.sin(2*np.pi*mesh.z/ .24)[None],(mesh.nr+1,mesh.nz))
    kernel=np.r_[radial.ravel(),np.zeros((mesh.nr-1)*mesh.nz),axial.ravel()]
    np.testing.assert_array_equal(difference_curl(g,kernel),0.)
    e=np.random.default_rng(382).normal(size=len(kernel))
    np.testing.assert_allclose(difference_curl(g,e),g.faraday@e,rtol=3e-14,atol=1e-10)


def test_kernel_separation_preserves_weak_physical_induction_without_cutoff():
    mesh=MHDMesh(.035,.12,5,10);model=CurrentThermal(mesh,np.full((6,10),1e18),39.948*atomic_mass);g=model.geometry
    rng=np.random.default_rng(241);signal=rng.normal(size=model.mass.shape[0]);signal/=np.linalg.norm(signal)
    radial=np.broadcast_to((1+mesh.r**2)[:,None],(mesh.nr,mesh.nz));axial=np.broadcast_to(np.cos(2*np.pi*mesh.z/.24)[None],(mesh.nr+1,mesh.nz))
    kernel=np.r_[radial.ravel(),np.zeros((mesh.nr-1)*mesh.nz),axial.ravel()]
    # Independently assemble the original unscaled three-block equations.
    C=g.faraday;H=model.hessian;M=model.mass;I=model.inertia;solver=InertialInduction(model)
    active=solver.active;Ca=C[active];Ha=H[active][:,active]
    matrix=bmat([[Ha,-Ha@Ca,None],[None,M,-I],[Ca.T@Ha,None,M]],format='csc')
    f=M@signal;right=np.r_[np.zeros(len(active)),f,np.zeros(M.shape[0])]
    direct=splu(matrix).solve(right)[:len(active)]
    for amplitude in (1.,1e-4,1e-8):
        actual=KernelInduction(model,tolerance=1e-10).solve(M@(kernel+amplitude*signal))
        error=actual.magnetic[active]/amplitude-direct
        relative=np.sqrt(error@(Ha@error)/(direct@(Ha@direct)))
        assert relative<max(3e-9,5e-14/amplitude),relative
        assert actual.residual<=1e-10
