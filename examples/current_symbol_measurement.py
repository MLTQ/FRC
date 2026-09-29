"""Reproducible combined Hall/pressure toroidal continuum-symbol evidence."""
from dataclasses import asdict
import numpy as np
from scipy.constants import electron_mass as ME, elementary_charge as E, mu_0 as MU0
from .electron_reference import ElectronReference
from .toroidal_thermal_symbol import coefficients,matrices,spectrum
from .toroidal_symbol_difference import direct_matrix


def measure():
    ref=ElectronReference();best=None
    radial=np.linspace(.001,.034,60);axial=np.linspace(-.12,.12,121);angle=np.linspace(0,np.pi,361)
    for r in radial:
        z,a=np.meshgrid(axial,angle,indexing='ij')
        _,_,_,aa,d,c,h,_=coefficients(ref,r,z,a)
        discriminant=(h-d)**2+4*aa*c;index=np.unravel_index(np.argmin(discriminant),discriminant.shape)
        value=float(discriminant[index])
        if best is None or value<best['discriminant_m2_s2']:
            best=dict(r_m=float(r),z_m=float(z[index]),angle_rad=float(a[index]),discriminant_m2_s2=value,
                      frozen_growth_per_wavenumber_m_s=float(.5*np.sqrt(max(0,-value))))
    radial_m=np.linspace(.002,.033,20);axial_m=np.linspace(-.12,.12,41);angle_m=np.linspace(0,np.pi,73)
    r,z,a=np.meshgrid(radial_m,axial_m,angle_m,indexing='ij');result=spectrum(ref,r,z,a)
    growth=np.max(abs(result['moving_eigenvalues'].imag),axis=-1);index=np.unravel_index(np.argmax(growth),growth.shape)
    location=tuple(float(v[index]) for v in (r,z,a));p,j,m,scale=matrices(ref,*location)
    physical=m*scale[None,:]/scale[:,None]
    derivatives=[]
    for spacing in (4e-6,2e-6,1e-6):
        direct=direct_matrix(ref,*location,scale,spacing=spacing)
        derivatives.append(dict(spacing_m=spacing,relative_error=float(np.linalg.norm(direct-physical)/np.linalg.norm(physical))))
    eig=result['moving_eigenvalues'][index]
    return dict(reference=asdict(ref),joint_grid=dict(radial_m=radial.tolist(),axial_m=axial.tolist(),angle_rad=angle.tolist()),
        joint_minimum=best,moving_grid=dict(radial_m=radial_m.tolist(),axial_m=axial_m.tolist(),angle_rad=angle_m.tolist()),
        moving_maximum=dict(r_m=location[0],z_m=location[1],angle_rad=location[2],growth_per_wavenumber_m_s=float(growth[index]),
            eigenvalue_real_m_s=eig.real.tolist(),eigenvalue_imaginary_m_s=eig.imag.tolist(),
            physical_matrix=m.tolist(),similarity_scales=scale.tolist(),normalized_matrix=physical.tolist()),
        nonreal_sample_count=int(np.sum(growth>1e-6)),total_samples=int(growth.size),
        independent_physical_derivatives=derivatives,physical_electron_skin_depth_m=float(np.sqrt(ME/(MU0*ref.density*E*E))),
        interpretation='Local frozen first-derivative toroidal sector of the stated massless isotropic closure; not a global well-posedness theorem or an attribution of the earlier axis failure.')
