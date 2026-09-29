"""Parent differential equations for unforced finite-mass compressible modes."""
from dataclasses import replace
import numpy as np
from scipy.constants import elementary_charge as E,electron_mass as ME,atomic_mass,mu_0 as MU0
from examples.finite_mass_wave_reference import FiniteMassWave


def wave(**kwargs):
    return FiniteMassWave(.035,.24,1e18,39.948*atomic_mass,.005,6e18*E,amplitude=1.,gamma_i=1.4,**kwargs)


def test_finite_mass_branches_satisfy_original_species_and_maxwell_equations():
    r,z=np.meshgrid(np.linspace(.002,.033,8),np.linspace(-.11,.1,9),indexing='ij')
    for ion,mass in ((39.948*atomic_mass,ME),(39.948*atomic_mass,0.),(ME,ME)):
        for branch in (0,1,2):
            ref=replace(wave(electron_mass=mass,branch=branch),ion_mass=ion);a=ref.fields(r,z,complex_output=True)
            h=1e-4*ref.radius
            def derivative(axis):
                fields=[ref.fields(r+i*h if axis==0 else r,z+i*h if axis==1 else z,complex_output=True) for i in (-2,-1,1,2)]
                return {k:(fields[0][k]-8*fields[1][k]+8*fields[2][k]-fields[3][k])/(12*h) for k in a}
            dr,dz,_=ref.derivatives(r,z);bias=np.array([0.,0.,ref.bias])
            norm=np.linalg.norm
            # Check the analytic derivative on its component scale, not the
            # tiny curl of the slow mode's large potential electric field.
            for actual,exact in ((derivative(0),dr),(derivative(1),dz)):
                for key in ('velocity','magnetic','electric'):
                    assert norm(actual[key]-exact[key])/norm(exact[key])<2e-9
            def curl(key):return np.stack((-dz[key][...,1],dz[key][...,0]-dr[key][...,2],dr[key][...,1]+a[key][...,1]/r),axis=-1)
            assert norm(curl('electric')-1j*ref.omega*a['magnetic'])/norm(ref.omega*a['magnetic'])<2e-8
            assert norm(curl('magnetic')/MU0-a['current'])/norm(a['current'])<2e-8
            n1=a['density']/(ref.ion_mass+ref.electron_mass)
            for name,ms,charge,gamma,p0 in (('ion',ref.ion_mass,E,ref.gamma_i,ref.ion_pressure),('electron',mass,-E,ref.gamma,ref.electron_pressure)):
                v=a[name+'_velocity'];gradient=a['pressure_gradient']*(ref.pressure_fraction if name=='electron' else 1-ref.pressure_fraction)
                inertial=-1j*ref.omega*ms*ref.number_density*v
                force=charge*ref.number_density*(a['electric']+np.cross(v,bias))-gradient
                scale=norm(charge*ref.number_density*a['electric'])+norm(gradient)+norm(inertial)
                assert norm(inertial-force)/scale<2e-11,(branch,name)
                div=dr[name+'_velocity'][...,0]+v[...,0]/r+dz[name+'_velocity'][...,2]
                derivative_scale=norm(dr[name+'_velocity'][...,0])+norm(v[...,0]/r)+norm(dz[name+'_velocity'][...,2])
                # Pair-plasma shear modes have exactly zero density/pressure
                # perturbation; normalize by original differential terms.
                assert norm(-1j*ref.omega*n1+ref.number_density*div)/(norm(ref.omega*n1)+ref.number_density*derivative_scale)<1e-10
                assert norm(-1j*ref.omega*a[name+'_pressure']+gamma*p0*div)/(norm(ref.omega*a[name+'_pressure'])+gamma*p0*derivative_scale)<1e-10
            lhs,rhs=ref.dispersion();x=ref.omega**2;v=ref.alfven_speed**2;s=ref.sound_speed**2;kk=ref.total_k**2
            scale=(ref.inertia*x+ref.k**2*v)*(ref.inertia*x*x+kk*(ref.inertia*s+v)*x+ref.k**2*kk*v*s)+abs(rhs)
            assert abs(lhs-rhs)/scale<2e-10
    # Equal physical masses cancel Hall without suppressing physical inertia.
    ref=replace(wave(),ion_mass=ME);assert ref.delta==0 and ref.inertia>1
    a=ref.fields(np.full(9,ref.radius),z[0],complex_output=True)
    assert np.linalg.norm(a['magnetic'][:,:2])<1e-14*np.linalg.norm(a['magnetic'])
