"""Representative uniform physical frequencies versus actual Gauss step sizes."""
import json
from pathlib import Path
import numpy as np
from .matched_mode_fixture import fixture
from .implicit_plasma_records import provenance,freeze,sha
from plasma_bench.modal_wall_coordinates import ModalWallCoordinates
from plasma_bench.modal_wall_inverse import ModalWallInverse


def calculate(spec,steps=(30,60,120),duration=1e-4):
    _,_,base,system,bias,_,_=fixture(spec)
    coordinates=ModalWallCoordinates(system,base,bias);inverse=ModalWallInverse(coordinates)
    spectra=[]
    for fluid,fb,bf,magnetic in inverse.operators:
        matrix=np.block([[fluid.toarray(),fb],[bf,magnetic]])
        spectra.append(np.linalg.eigvals(matrix))
    spectra=np.array(spectra);frequency=abs(spectra.imag)
    nonzero=frequency>frequency.max()*1e-10;omega=frequency[nonzero]
    records=[]
    for count in steps:
        dt=duration/count;x=dt*omega
        # Continuous branch from zero to 2*pi, not a principal-phase alias.
        phase=2*np.arctan2(x/2,1-x*x/12)
        records.append(dict(steps=count,dt=dt,maximum_omega_dt=float(x.max()),
            fraction_nonzero_modes_with_omega_dt_above_one=float(np.mean(x>1)),
            maximum_unwrapped_frequency_relative_error=float(np.max(abs(phase/x-1))),
            maximum_uniform_mode_endpoint_error=float(np.max(abs(np.exp(1j*count*phase)-np.exp(1j*duration*omega))))))
    return dict(spec=spec,duration=duration,maximum_frequency_radians_per_second=float(omega.max()),
        maximum_real_part_over_frequency=float(np.max(abs(spectra.real))/omega.max()),
        nonzero_mode_count=int(nonzero.sum()),time_controls=records),dict(spectra=spectra,background=coordinates.background())


def run(destination,specs):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    identity=provenance(['examples.matched_mode_phase_audit']);destination.mkdir(parents=True);freeze(destination,identity)
    cases=[]
    for spec in specs:
        record,arrays=calculate(spec);raw=destination/f'case_{len(cases):02}.npz';np.savez_compressed(raw,**arrays)
        record.update(arrays_file=raw.name,arrays_sha256=sha(raw));cases.append(record)
        print('representative uniform phase',spec,record['maximum_frequency_radians_per_second'],record['time_controls'],flush=True)
    if provenance(['examples.matched_mode_phase_audit'])!=identity:raise RuntimeError('Source changed')
    result=dict(terminal=True,source_sha256=identity['source_sha256'],cases=cases,
        limitations=['Uniform representative operator, not the actual nonuniform spectrum or a bound',
            'No physical amplitude projection or diagnosis of which transient was excited',
            'Phase underresolution explains a risk; full actual response convergence remains required'])
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
