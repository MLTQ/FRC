"""Independent annular moments and smooth nearly solenoidal flux checks."""
import numpy as np
from scipy.constants import atomic_mass,elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.annular_polynomial import radial_matrix,axial_samples,normal_momentum_flux
from examples.finite_mass_wave_reference import FiniteMassWave
from examples.finite_mass_manufactured import cell_average


def test_cubic_annular_average_interpolation_is_exact_through_axis_and_wall():
    mesh=MHDMesh(.035,.12,24,16);r=np.linspace(0,mesh.radius,71)
    for power in (0,1,2):
        averages=np.zeros(mesh.nr);exact=np.zeros_like(r)
        for degree,coefficient in enumerate((1.,-.7,.4,-.2)):
            exponent=power+2*degree
            averages+=coefficient*(mesh.rf[1:]**(exponent+2)-mesh.rf[:-1]**(exponent+2))/((exponent+2)*mesh.radial_measure*mesh.radius**exponent)
            exact+=coefficient*(r/mesh.radius)**exponent
        np.testing.assert_allclose(radial_matrix(mesh,r,power=power)@averages,exact,rtol=3e-13,atol=3e-14)
    # Fourth-order periodic sine reconstruction from analytic cell averages.
    errors=[]
    for nz in (16,32,64):
        h=2*np.pi/nz;faces=np.arange(nz+1)*h
        averages=(np.cos(faces[:-1])-np.cos(faces[1:]))/h
        errors.append(np.linalg.norm(axial_samples(averages[None],.23)[0]-np.sin(faces[:-1]+.23*h)))
    assert errors[1]<errors[0]/10 and errors[2]<errors[1]/10


def test_nearly_solenoidal_mode_uses_true_normal_momentum_averages():
    ref=FiniteMassWave(.035,.24,1e18,39.948*atomic_mass,.005,6e18*E,gamma_i=1.4,amplitude=1.,branch=0)
    errors=[]
    for nr in (8,16,32):
        mesh=MHDMesh(ref.radius,ref.length/2,nr,2*nr)
        def linear_momentum(r,z):
            f=ref.fields(r,z);v=f['velocity']
            return np.stack((np.zeros_like(r),ref.density*v[...,0],ref.density*v[...,2],np.zeros_like(r),ref.density*r*v[...,1]),axis=-1)
        averages=cell_average(mesh,linear_momentum);radial,axial=normal_momentum_flux(mesh,averages)
        numerical=-(np.diff(radial,axis=0)+np.diff(axial,axis=1))/mesh.volumes
        exact=cell_average(mesh,lambda r,z:ref.derivatives(r,z)[2]['density'].real)
        errors.append(np.sqrt(np.sum(mesh.volumes*(numerical-exact)**2)/np.sum(mesh.volumes*exact**2)))
    assert errors[1]<errors[0]/7 and errors[2]<errors[1]/7,errors
