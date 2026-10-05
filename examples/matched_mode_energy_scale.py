"""Independent small-signal normalization of actual trajectory energy drift."""
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.regular_current import RegularCurrent
from .matched_mode_fixture import fixture
from .matched_mode_comparison import load
from .current_transport_sampling import TransportSampler
from .matched_mode_rates import fingerprint
from .perturbation_energy_reference import densities
from .implicit_plasma_records import provenance,freeze,sha


def measure(spec,epsilon,order):
    mode,d,base,system,bias,du,db=fixture(spec)
    s=TransportSampler(d,system.geometry.mesh,order)
    n,_,_,dn,dv=s.fluid(base,du)
    entropy=np.stack([s.scalar(base[...,k])/n for k in (4,5)],axis=-1)
    ds=np.stack([(s.scalar(du[...,k])-entropy[...,k-4]*dn)/n for k in (4,5)],axis=-1)
    physical=epsilon*np.concatenate((dn[...,None],dv,ds,s.magnetic(db)),axis=-1)
    model=RegularCurrent(system.mapping,base[...,0],system.ion_mass)
    current=epsilon*s.electric(model.current(db))
    pressure=E*n[...,None]*np.exp(entropy)*(n[...,None]/d.thermal.n0)**(d.thermal.gammas-1)
    actual=densities(n,pressure,d.thermal.gammas,d.mass,model.reduced_mass,physical,current)
    reference=mode.fields(s.R,s.Z);background=mode.ref.fields(s.Z)
    reference_pressure=np.stack((background['electron_pressure'],background['ion_pressure']),axis=-1)
    expected=densities(background['number_density'],reference_pressure,d.thermal.gammas,d.mass,model.reduced_mass,
        epsilon*reference['physical'],epsilon*reference['current'])
    integrals={name:float(np.sum(s.weight*value)) for name,value in actual.items()}
    continuum={name:float(np.sum(s.weight*value)) for name,value in expected.items()}
    total=sum(integrals.values());target=sum(continuum.values())
    if not np.isfinite([total,target]).all() or min(total,target)<=0:raise ValueError('Positive finite physical signal scale required')
    return dict(order=order,quadratic_energy=total,continuum_quadratic_energy=target,
        component_energies=integrals,continuum_component_energies=continuum,
        relative_representation_difference=(total-target)/target,
        physical_fingerprints={name:fingerprint(value) for name,value in actual.items()},
        continuum_fingerprints={name:fingerprint(value) for name,value in expected.items()}),dict(base_state=base,fluid_direction=du,magnetic_direction=db)


def run(audits,destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    identity=provenance(['examples.matched_mode_energy_scale']);destination.mkdir(parents=True);freeze(destination,identity)
    cases=[]
    for audit_directory in audits:
        audit,arrays=load(audit_directory);quadrature=[]
        for original in audit['quadrature']:
            value,raw=measure(audit['spec'],audit['epsilon'],original['order'])
            value['signed_energy_drift_over_initial_quadratic']={name:record['energy_change']/value['quadratic_energy']
                for name,record in original['conservation'].items()}
            value['initial_quadratic_over_background']={name:value['quadratic_energy']/record['initial_energy']
                for name,record in original['conservation'].items()}
            quadrature.append(value)
        path=destination/f'case_{len(cases):02}.npz';np.savez_compressed(path,**raw)
        cases.append(dict(audit_directory=str(audit_directory),input_summary_sha256=sha(Path(audit_directory)/'summary.json'),
            spec=audit['spec'],epsilon=audit['epsilon'],steps=audit['steps'],duration=audit['duration'],quadrature=quadrature,
            arrays_file=path.name,arrays_sha256=sha(path)))
        print('signal-energy normalization',audit_directory,quadrature[-1]['signed_energy_drift_over_initial_quadratic'],flush=True)
    if provenance(['examples.matched_mode_energy_scale'])!=identity:raise RuntimeError('Source changed')
    result=dict(terminal=True,source_sha256=identity['source_sha256'],cases=cases,
        definition='Initial physical energy Hessian on intended conservative perturbation; scale only, not an invariant or repair')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
