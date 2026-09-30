import numpy as np
import pytest
from scipy.constants import atomic_mass,electron_mass
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.compatible_density import CompatibleDensity
from plasma_bench.conservative_kinematics import ConservativeKinematics
from plasma_bench.log_fluid import LogFluid
from plasma_bench.compatible_particle_flux import CompatibleParticleFlux
from examples.log_fluid_reference import LogFluidReference,average,initial


def test_independent_nonlinear_characteristics_and_rotating_euler():
    r=np.linspace(.002,.032,8)[:,None];z=np.linspace(-.113,.117,11)[None,:]
    for kind,direction in [('contact',1),('contact',-1),('acoustic',1),('acoustic',-1),('rotation',1)]:
        ref=LogFluidReference(kind=kind,direction=direction,gamma_i=5/3 if kind=='acoustic' else 1.4)
        f=ref.fields(r,z,.8*ref.duration);dr,dz,dt=ref.derivatives(r,z,.8*ref.duration)
        v=f['velocity'];div=dr['velocity'][...,0]+v[...,0]/r+dz['velocity'][...,2]
        continuity=dt['n']+v[...,0]*dr['n']+v[...,2]*dz['n']+f['n']*div
        assert np.max(abs(continuity))/(ref.number_density/ref.duration)<2e-13
        acceleration=dt['velocity']+v[...,0,None]*dr['velocity']+v[...,2,None]*dz['velocity']
        acceleration[...,0]-=v[...,1]**2/r
        acceleration[...,1]+=v[...,0]*v[...,1]/r
        pressure_gradient=np.stack((dr['pe']+dr['pi'],np.zeros_like(div),dz['pe']+dz['pi']),axis=-1)
        residual=ref.mass*f['n'][...,None]*acceleration+pressure_gradient
        assert np.max(abs(residual))/(ref.pressure/ref.radius)<3e-13
        for key,gamma in [('pe',ref.gamma_e),('pi',ref.gamma_i)]:
            residual=dt[key]+v[...,0]*dr[key]+v[...,2]*dz[key]+gamma*f[key]*div
            assert np.max(abs(residual))/(ref.pressure/ref.duration)<3e-13


def test_true_kinetic_gradients_and_regular_curl_converge():
    errors=[];axis=[]
    for nr in (8,16,32):
        mesh=MHDMesh(.035,.12,nr,2*nr);geo=CompatibleGeometry(mesh,degree=5,quadrature_order=12)
        M=39.948*atomic_mass+electron_mass;k=np.pi/.12
        def fields(r,z):
            y=(r/.035)**2;n=1e18*(1+.3*y)*(1+.15*np.cos(k*z))
            v=np.stack((200*r/.035*(1-y)*(1+.2*np.cos(k*z)),300*r/.035*(1+.1*y)*np.sin(k*z+.2),400*(1+.2*y)*np.cos(k*z)),axis=-1)
            m=np.stack((n*v[...,0],n*v[...,2],n*r*v[...,1]),axis=-1)
            return n,v,m
        density=CompatibleDensity(geo,average(mesh,lambda r,z:fields(r,z)[0]))
        momentum=average(mesh,lambda r,z:fields(r,z)[2]);kin=ConservativeKinematics(geo,M)
        f=kin.fields(density,momentum);grad=kin.gradients(density,momentum)
        n,v,_=fields(geo.r,geo.z);h=1e-28
        vr=np.imag(fields(geo.r+1j*h,geo.z)[1])/h;vz=np.imag(fields(geo.r,geo.z+1j*h)[1])/h
        exact=np.stack((-vz[...,1],vz[...,0]-vr[...,2],vr[...,1]+v[...,1]/geo.r),axis=-1)
        err=np.sum((grad['curl']-exact)**2,axis=-1);norm=np.sum(exact**2,axis=-1)
        errors.append(np.sqrt(np.sum(geo.weight*err)/np.sum(geo.weight*norm)))
        axis.append(np.sqrt(np.sum(geo.weight[0]*err[0])/np.sum(geo.weight[0]*norm[0])))
        nt=average(mesh,lambda r,z:fields(r,z)[0]*.1*np.cos(k*z))
        mt=momentum*.3;rate=kin.energy_rate(density,momentum,nt,mt)
        derivatives=kin.energy_gradients(density,momentum)
        actual=np.sum(mesh.volumes*rate['total']);pull=np.sum(derivatives['density']*nt)+np.sum(derivatives['momentum']*mt)
        assert abs(actual-pull)/abs(actual)<2e-13
        expected=.5*M*n*np.sum(v*v,axis=-1)
        assert abs(np.sum(geo.weight*(f['kinetic']-expected)))/np.sum(geo.weight*expected)<2e-5
    assert errors[0]/errors[1]>20 and errors[1]/errors[2]>20
    assert axis[0]/axis[1]>20 and axis[1]/axis[2]>20


def make(nr,kind='contact'):
    ref=LogFluidReference(kind=kind,gamma_i=5/3 if kind=='acoustic' else 1.4)
    geo=CompatibleGeometry(MHDMesh(ref.radius,ref.half_length,nr,2*nr),degree=5,quadrature_order=12)
    driver=LogFluid(geo,ref.ion_mass,gamma_e=ref.gamma_e,gamma_i=ref.gamma_i)
    return ref,driver,initial(driver,ref)


def test_shared_continuity_and_conservative_momentum_energy_fluxes():
    for kind in ('contact','acoustic','rotation'):
        ref,driver,state=make(8,kind);rate,d=driver.derivative(state)
        v=driver.mesh.volumes
        assert d['continuity_face_agreement']<3e-13
        model,species,*_=driver.evaluate(state)
        particle=CompatibleParticleFlux(model.density_map)
        nf=particle.integrated(particle.from_momentum(state.momentum[...,:2]))
        mechanical=driver.transport.rate(model,driver.kinematics,state.momentum,species,(state.electron,state.ion))
        fluxes=[nf]+[(mechanical['radial_momentum_flux'][...,k]/ref.mass,
                      mechanical['axial_momentum_flux'][...,k]/ref.mass) for k in (1,2)]
        fluxes.append((mechanical['radial_energy_flux'],mechanical['axial_energy_flux']))
        for q,faces in zip((rate.density,rate.momentum[...,1],rate.momentum[...,2],rate.energy_account),fluxes):
            # Backward error uses the uncancelled flux, including stationary cases.
            scale=sum(np.sum(abs(f)) for f in faces)
            assert abs(np.sum(v*q))<1e-12*max(scale,1e-30)
        assert np.isfinite(d['local_energy_defect']).all()


def test_unforced_nonlinear_wave_refines_and_failure_is_atomic():
    errors=[]
    for nr in (8,16):
        ref,driver,state=make(nr,'acoustic');end=.05*ref.duration;steps=12
        for i in range(steps):state,_=driver.step(state,i*end/steps,end/steps)
        n=driver.evaluate(state)[0].density;g=driver.geometry;exact=ref.fields(g.r,g.z,end)['n']
        errors.append(np.sqrt(np.sum(g.weight*(n-exact)**2)/np.sum(g.weight*(exact-ref.number_density)**2)))
    assert errors[0]/errors[1]>15
    ref,driver,state=make(8,'acoustic');before=[v.copy() for v in (state.density,state.momentum,state.electron,state.ion,state.energy_account)]
    with pytest.raises((ValueError,RuntimeError)):
        driver.step(state,0.,100*ref.duration)
    for a,b in zip(before,(state.density,state.momentum,state.electron,state.ion,state.energy_account)):assert np.array_equal(a,b)
