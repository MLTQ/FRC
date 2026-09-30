import numpy as np
from plasma_bench.radial_pressure_fluid import RadialPressureFluid
from plasma_bench.compatible_particle_flux import CompatibleParticleFlux,ParticleFaces
from plasma_bench.compatible_bulk_thermal import CompatibleBulkThermal
from examples.log_fluid_study import initialize
from examples.log_fluid_reference import initial
from examples.radial_work_reference import RadialWorkReference


def test_radial_advection_and_pressure_power_pair_with_true_kinetic_energy():
    ref,base,_=initialize(dict(kind='contact',nr=8));ref=RadialWorkReference()
    d=RadialPressureFluid(base.geometry,ref.ion_mass,gamma_e=ref.gamma_e,gamma_i=ref.gamma_i);state=initial(d,ref)
    model,species,*_=d.evaluate(state);flux=CompatibleParticleFlux(model.density_map)
    faces=flux.from_momentum(state.momentum[...,:2]);radial_faces=ParticleFaces(faces.radial,np.zeros_like(faces.axial))
    rates=[CompatibleBulkThermal(t,flux).rate(s,radial_faces) for t,s in zip(species,(state.electron,state.ion))]
    d.derivative(state)
    work=sum(t.face_work.radial for t in rates)+d.radial_kinetic_work(model,flux,state)
    source,residual=d.radial_pressure.solve(model,flux,work)
    actual=d.transport.rate(model,d.kinematics,state.momentum,species,(state.electron,state.ion))
    pressure=d.transport.rate(model,d.kinematics,np.zeros_like(state.momentum),species,(state.electron,state.ion))
    mr=-(np.diff(actual['radial_momentum_flux']-pressure['radial_momentum_flux'],axis=0)/d.mesh.volumes[...,None])/d.mass
    mr[...,0]+=source;nt=flux.density_rate(radial_faces)
    kinetic=d.kinematics.energy_rate(model.density_map,state.momentum,nt,mr)['total']
    heat=[t.energy_rate(s,k.entropy,nt)['total'] for t,s,k in zip(species,(state.electron,state.ion),rates)]
    scale=np.sum(d.mesh.volumes*(abs(kinetic)+sum(abs(q) for q in heat)))
    assert abs(np.sum(d.mesh.volumes*(kinetic+sum(heat))))/scale<5e-12
    assert residual<2e-12


def test_two_dimensional_parent_rates_and_rotating_balance_refine():
    errors=[];rotation=[]
    for nr in (8,16,32):
        ref,b,s=initialize(dict(kind='rotation',nr=nr));d=RadialPressureFluid(b.geometry,ref.ion_mass,gamma_e=ref.gamma_e,gamma_i=ref.gamma_i)
        rate,_=d.derivative(s);model=d.evaluate(s)[0];g=d.geometry
        acceleration=d.kinematics.sampled_momentum(model.density_map,rate.momentum)[...,0]/model.density
        rotation.append(np.sqrt(np.sum(g.weight*acceleration**2)/d.mesh.volumes.sum()))
        ref=RadialWorkReference();state=initial(d,ref);model=d.evaluate(state)[0];rate,_=d.derivative(state)
        n=model.density;nt=model.density_map.sample(rate.density);v=d.kinematics.fields(model.density_map,state.momentum)['velocity']
        mt=d.kinematics.sampled_momentum(model.density_map,rate.momentum);vt=(mt-v*nt[...,None])/n[...,None]
        target=ref.rates(g.r,g.z);err=np.sqrt(np.sum(g.weight[...,None]*(vt-target['velocity'])**2)/np.sum(g.weight[...,None]*target['velocity']**2))
        errors.append(err)
    assert all(a/b>10 for a,b in zip(rotation[:-1],rotation[1:]))
    assert errors[0]/errors[1]>3 and errors[1]/errors[2]>3
