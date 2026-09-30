"""Physical plane-wave transport, all-mode growth and annular face accuracy."""
import numpy as np
from scipy.linalg import expm
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.normal_momentum_faces import axial_faces,radial_faces
from plasma_bench.particle_flux_basis import difference
from examples.finite_mass_manufactured import cell_average


def test_periodic_uniform_translation_has_no_semidiscrete_or_rk4_mode_growth():
    for degree in (3,5):
        mesh=MHDMesh(.035,.12,8,32)
        for speed in (-2000.,2000.):
            operator=-speed*(difference(mesh.zf,True)@axial_faces(mesh,degree)).toarray()
            np.testing.assert_allclose(operator+operator.T,0,atol=2e-9)
            eigenvalues=np.linalg.eigvals(operator);dt=.2*mesh.dz/abs(speed);z=dt*eigenvalues
            amplification=1+z+z*z/2+z**3/6+z**4/24
            assert max(abs(amplification))<=1+2e-13
            values=np.random.default_rng(21).normal(size=mesh.nz)
            evolved=expm(40*(2*mesh.half_length/abs(speed))*operator)@values
            np.testing.assert_allclose(np.linalg.norm(evolved),np.linalg.norm(values),rtol=3e-11)


def test_cell_average_plane_wave_advection_and_regular_radial_moments():
    for degree in (3,5):
        errors=[]
        for nz in (16,32):
            mesh=MHDMesh(.035,.12,8,nz);k=np.pi/mesh.half_length
            averages=np.diff(-np.cos(k*mesh.zf))/(k*mesh.dz)
            exact=-2000*np.diff(np.sin(k*mesh.zf))/mesh.dz
            result=-2000*difference(mesh.zf,True)@axial_faces(mesh,degree)@averages
            errors.append(np.linalg.norm(result-exact)/np.linalg.norm(exact))
        assert errors[1]<errors[0]/(12 if degree==3 else 45),errors
        for power in (0,1,2):
            mesh=MHDMesh(.035,.12,12,24)
            def physical(r,z):return (r/mesh.radius)**power*(.3+.4*(r/mesh.radius)**2-.2*(r/mesh.radius)**(2*degree))
            means=cell_average(mesh,physical,order=16)
            value=radial_faces(mesh,degree,power=power)@means
            target=physical(mesh.rf[:,None],mesh.z[None])
            np.testing.assert_allclose(value,np.broadcast_to(target,value.shape),rtol=2e-11,atol=2e-13)
