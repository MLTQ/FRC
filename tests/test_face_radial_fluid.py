from dataclasses import fields
import numpy as np
import pytest
from plasma_bench.compatible_particle_flux import CompatibleParticleFlux,ParticleFaces
from plasma_bench.compatible_bulk_thermal import CompatibleBulkThermal
from plasma_bench.face_radial_fluid import FaceRadialState
from examples.face_radial_reference import initialize
from examples.face_radial_norms import FaceRadialSampler
from examples.radial_work_reference import RadialWorkReference


def test_regular_polynomial_force_and_true_face_kinetic_work():
    ref,d,s=initialize(dict(nr=8));model=d.evaluate(s)[0];g=d.geometry
    # Nonuniform actual n is intentional: physical force is n*a, not a fitted mass.
    flux=CompatibleParticleFlux(model.density_map)
    acceleration=g.r/ref.radius*(1-(g.r/ref.radius)**2)**2
    values=np.zeros(g.r.shape+(3,));values[...,0]=-d.mass*acceleration
    work=flux.vector_adjoint(values).radial
    qt,residual=d.kinematics.solve_rate(model.density_map,np.zeros_like(g.r),work)
    actual=d.kinematics.sampled_momentum(model.density_map,qt,np.zeros_like(s.momentum))[...,0]/model.density
    # This contact density varies only in z, so n*a belongs to the declared space.
    assert np.sqrt(np.sum(g.weight*(actual-acceleration)**2)/np.sum(g.weight*acceleration**2))<2e-11
    assert np.sqrt(np.sum(g.weight[0]*(actual[0]-acceleration[0])**2)/np.sum(g.weight[0]*acceleration[0]**2))<2e-11
    assert residual<2e-12
    target=RadialWorkReference();_,d,s=initialize(dict(nr=8),target);model=d.evaluate(s)[0];g=d.geometry
    nref=target.fields(g.r,g.z)['n'];vref=target.fields(g.r,g.z)['velocity']
    k=d.kinematics.fields(model.density_map,s.radial,s.momentum)['kinetic']
    assert abs(np.sum(g.weight*(k-.5*d.mass*nref*np.sum(vref*vref,axis=-1))))/np.sum(g.weight*k)<2e-4


def test_diagnostic_radial_work_identity_and_actual_cell_momenta_conserve():
    ref,d,s=initialize(dict(nr=8),RadialWorkReference());model,species,_,kinetic,_=d.evaluate(s)
    flux=CompatibleParticleFlux(model.density_map);f=d.particle_faces(flux,s)
    faces=ParticleFaces(f.radial,np.zeros_like(f.axial));nt=flux.density_rate(faces)
    thermal=[CompatibleBulkThermal(t,flux).rate(v,faces) for t,v in zip(species,(s.electron,s.ion))]
    mechanical=d.transport.rate(model,s.radial,s.momentum,species,(s.electron,s.ion))
    grad=d.kinematics.gradients(model.density_map,s.radial,s.momentum);particle=flux.sample(faces)
    advection=-particle[...,0]*grad['radial'][...,0]+kinetic['velocity'][...,0]*model.density_map.sample(nt)
    work=sum(t.face_work.radial for t in thermal)+d.radial_kinetic_work(model,flux,s,mechanical,kinetic,grad)
    qt,_=d.kinematics.solve_rate(model.density_map,advection,work)
    mt=-np.diff(mechanical['radial_momentum_flux'],axis=0)/(d.mass*d.mesh.volumes[...,None])
    k=d.kinematics.energy_rate(model.density_map,s.radial,s.momentum,nt,qt,mt)['total']
    heat=[t.energy_rate(v,r.entropy,nt)['total'] for t,v,r in zip(species,(s.electron,s.ion),thermal)]
    scale=np.sum(d.mesh.volumes*(abs(k)+sum(abs(v) for v in heat)))
    assert abs(np.sum(d.mesh.volumes*(k+sum(heat))))/scale<5e-12
    rate=FaceRadialState(nt,qt,mt,thermal[0].entropy,thermal[1].entropy,np.zeros_like(nt))
    _,independent=FaceRadialSampler(d.mesh,order=16).rates(s,rate,ref,0.)
    assert abs(np.sum(d.mesh.volumes*independent))/scale<5e-12
    full,diagnostic=d.derivative(s)
    assert diagnostic['continuity_face_agreement']<2e-12
    for component in (0,1):
        val=full.momentum[...,component]
        assert abs(np.sum(d.mesh.volumes*val))/np.sum(d.mesh.volumes*abs(val))<2e-12


def test_original_two_dimensional_rates_and_axis_force_refine():
    errors=[];axis=[];curl=[]
    for nr in (8,16,32):
        ref,d,s=initialize(dict(nr=nr),RadialWorkReference());rate,_=d.derivative(s);model=d.evaluate(s)[0];g=d.geometry
        norms,_=FaceRadialSampler(d.mesh).rates(s,rate,ref,0.)
        errors.append(norms['velocity']['relative_L2']);axis.append(norms['velocity_components']['radial']['axis_relative_L2'])
        h=1e-28;v=ref.fields(g.r,g.z)['velocity']
        dr=np.imag(ref.fields(g.r+1j*h,g.z)['velocity'])/h;dz=np.imag(ref.fields(g.r,g.z+1j*h)['velocity'])/h
        exact=np.stack((-dz[...,1],dz[...,0]-dr[...,2],dr[...,1]+v[...,1]/g.r),axis=-1)
        actual=d.kinematics.gradients(model.density_map,s.radial,s.momentum)['curl']
        curl.append(np.sqrt(np.sum(g.weight[...,None]*(actual-exact)**2)/np.sum(g.weight[...,None]*exact**2)))
    assert errors[-1]<1e-3 and axis[-1]<1e-3
    assert all(a/b>8 for a,b in zip(errors[:-1],errors[1:]))
    assert all(a/b>8 for a,b in zip(axis[:-1],axis[1:]))
    assert all(a/b>15 for a,b in zip(curl[:-1],curl[1:]))


def test_unforced_simple_wave_refines_and_failed_step_is_atomic():
    errors=[]
    for nr in (8,16):
        ref,d,s=initialize(dict(nr=nr,kind='acoustic'));end=.05*ref.duration;steps=12
        for i in range(steps):s,_=d.step(s,i*end/steps,end/steps)
        n=d.evaluate(s)[0].density;g=d.geometry;exact=ref.fields(g.r,g.z,end)['n']
        errors.append(np.sqrt(np.sum(g.weight*(n-exact)**2)/np.sum(g.weight*(exact-ref.number_density)**2)))
    assert errors[0]/errors[1]>15
    ref,d,s=initialize(dict(nr=8,kind='acoustic'));before={f.name:getattr(s,f.name).copy() for f in fields(s)}
    with pytest.raises((ValueError,RuntimeError)):d.step(s,0.,100*ref.duration)
    for key,value in before.items():assert np.array_equal(getattr(s,key),value)


def test_physical_pressure_force_pairs_only_thermal_work():
    ref,d,s=initialize(dict(nr=8),RadialWorkReference());model,species,_,_,_=d.evaluate(s)
    flux=CompatibleParticleFlux(model.density_map);all_faces=d.particle_faces(flux,s)
    faces=ParticleFaces(all_faces.radial,np.zeros_like(all_faces.axial));nt=flux.density_rate(faces)
    thermal=[CompatibleBulkThermal(t,flux).rate(v,faces) for t,v in zip(species,(s.electron,s.ion))]
    qt,_=d.kinematics.solve_rate(model.density_map,np.zeros_like(d.geometry.r),sum(t.face_work.radial for t in thermal))
    kt=d.kinematics.energy_rate(model.density_map,s.radial,s.momentum,np.zeros_like(nt),qt,np.zeros_like(s.momentum))['total']
    ht=sum(t.energy_rate(v,k.entropy,nt)['total'] for t,v,k in zip(species,(s.electron,s.ion),thermal))
    scale=np.sum(d.mesh.volumes*(abs(kt)+abs(ht)))
    assert abs(np.sum(d.mesh.volumes*(kt+ht)))/scale<5e-12


def test_interrupted_nonlinear_wave_resumes_identical_accepted_history():
    from examples.face_radial_study import integrate
    from examples.face_radial_records import save,sha
    import json
    import tempfile
    from pathlib import Path
    spec=dict(kind='acoustic',nr=8,duration_factor=.02);saved=[]
    class Interrupted(Exception):pass
    def checkpoint(record,arrays):
        saved.append((record,{k:np.asarray(v).copy() for k,v in arrays.items()}));raise Interrupted
    with pytest.raises(Interrupted):integrate(spec,checkpoint)
    with tempfile.TemporaryDirectory() as temp:
        path=Path(temp);save(path,'fluid_00',*saved[0]);record=json.loads((path/'fluid_00_checkpoint.json').read_text())
        assert sha(path/record['arrays_file'])==record['arrays_sha256']
        with np.load(path/record['arrays_file']) as data:resume=(record,{k:data[k].copy() for k in data.files})
        completed,restarted=integrate(spec,resume=resume)
    reference,uninterrupted=integrate(spec)
    assert completed['status']==reference['status']=='complete'
    assert completed['accepted_steps']==reference['accepted_steps']
    assert set(restarted)==set(uninterrupted)
    assert all(np.array_equal(restarted[k],v) for k,v in uninterrupted.items())
