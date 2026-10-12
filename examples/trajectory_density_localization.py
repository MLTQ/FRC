"""Radial, axial and interface localization of independently audited density."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import numpy as np
from scipy.special import j0
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.mhd_mesh import MHDMesh
from .continuum_contact_mode import ContinuumContactMode
from .current_transport_sampling import TransportSampler
from .trajectory_density_prefix import fingerprint
from .implicit_plasma_records import provenance, freeze, sha


def calculate(case, observed):
    directory = Path(case['directory']); record = case['summary_snapshot']; spec = record['spec']
    steps = case['prefix_steps']
    if steps != observed['prefix_steps'] or str(directory) != observed['directory']:
        raise ValueError('Matching previously audited prefix required')
    inputs = observed['input_raw_sha256']
    if any(sha(directory/p) != h for p,h in inputs.items()):
        raise ValueError('Original accepted state/stage bytes required')
    source = json.loads((directory/'provenance.json').read_text())
    if source['source_sha256'] != observed['physical_source_sha256']:
        raise ValueError('Original physical source required')
    if any(sha(p) != h for p,h in source['source_files'].items()):
        raise ValueError('Physical source changed')
    mode = ContinuumContactMode(); ref = mode.ref
    if mode.input_sha256 != observed['reference_input_sha256']:
        raise ValueError('Original independent continuum mode required')
    fluid = RegularCoordinateFluid(ref.radius,ref.length,spec['elements'],spec['points'],
        ref.ion_mass+ref.electron_mass,degree=5,gammas=(ref.gamma_e,ref.gamma_i))
    mesh = MHDMesh(ref.radius,ref.length/2,spec['field_nr'],spec['field_nz'])
    if mesh.nr % fluid.elements: raise ValueError('Field cells must subdivide fluid elements')
    with np.load(directory/'initial.npz') as raw: initial = raw['state'][...,0].copy()
    selected = record['records'][:steps]
    names = ['initial.npz',*[r['arrays_file'] for r in selected]]
    times = np.array([0.,*[r['time'] for r in selected]])
    wave = np.rint(np.fft.fftfreq(mesh.nz)*mesh.nz).astype(int)
    interfaces = ref.radius*np.arange(1,fluid.elements)/fluid.elements
    arrays = dict(time=times,signed_cell_mode=wave,interface_radius=interfaces,
        fluid_radial_edges=np.linspace(0,ref.radius,fluid.elements+1)); quadrature = []
    for order in (16,20):
        original = next(q for q in observed['quadrature'] if q['order']==order)
        sample = TransportSampler(fluid,mesh,order)
        direction = j0(mode.q*sample.R)*mode.amplitudes(sample.Z)[...,0]
        intended = record['epsilon']*direction
        base = ref.fields(sample.Z)['number_density']
        initial_error = sample.scalar(initial)-base-intended
        if fingerprint(sample.weight)!=original['weight_sha256'] or fingerprint(direction)!=original['reference_sha256']:
            raise ValueError('Original physical integration and reference fingerprints required')
        axial_weight = np.tile(mesh.dz/2*np.polynomial.legendre.leggauss(order)[1],mesh.nz)
        surface_weight = 2*np.pi*interfaces[:,None]*axial_weight[None,:]
        trace_base_energy = np.sum(surface_weight*base*base,axis=1)
        traces=[]; radial=[]; spectra=[]; metrics=[]; checks=[]
        for index,(filename,when) in enumerate(zip(names,times)):
            with np.load(directory/filename) as raw: density = raw['state'][...,0].copy()
            delta = sample.scalar(density-initial)
            error = delta-np.expm1(mode.growth*when)*intended+initial_error
            hashes = dict(density_change=fingerprint(delta),reference_error=fingerprint(error))
            if hashes != original['physical_fingerprints'][index]:
                raise ValueError('Original complete physical density samples required')
            squared = error*error*sample.weight
            total = float(np.sum(squared))
            partition = np.sum(squared.reshape(fluid.elements,-1,len(sample.z)),axis=(1,2))
            transformed = np.fft.fft(error.reshape(len(sample.r),mesh.nz,order),axis=1)/mesh.nz
            spectrum = mesh.nz*np.sum(abs(transformed)**2*sample.weight[:,None,:order],axis=(0,2))
            defects = [abs(float(np.sum(value))-total)/total if total else 0. for value in (partition,spectrum)]
            if max(defects)>256*np.finfo(float).eps:
                raise ValueError('Complete radial partition or discrete Parseval identity failed')
            jump = density[:-1,-1,:]-density[1:,0,:]
            trace = (np.fft.fft(jump,axis=1)/fluid.nz)@sample.fourier.T
            if np.max(abs(trace.imag),initial=0.)>1e-10*max(np.max(abs(trace.real),initial=0.),np.finfo(float).tiny):
                raise ValueError('Nonreal independently reconstructed interface jump')
            trace = trace.real
            jump_energy = np.sum(surface_weight*trace*trace,axis=1)
            trace_relative = np.sqrt(jump_energy/trace_base_energy)
            radial.append(partition); spectra.append(spectrum); traces.append(jump_energy)
            metrics.append(dict(time=float(when),error_squared_integral=total,
                relative_L2_error=float(np.sqrt(total))/(np.exp(mode.growth*when)*original['initial_signal_L2']),
                largest_error_element=int(np.argmax(partition)),
                element_error_fractions=(partition/total if total else partition).tolist(),
                largest_signed_cell_mode=int(wave[np.argmax(spectrum)]),
                low_mode_fraction=float(np.sum(spectrum[abs(wave)<=4])/total) if total else 0.,
                upper_quarter_mode_fraction=float(np.sum(spectrum[abs(wave)>=3*mesh.nz/8])/total) if total else 0.,
                interface_jump_L2_over_background=trace_relative.tolist(),
                sampled_max_interface_jump_over_background=float(np.max(abs(trace/base),initial=0.)),
                partition_relative_defect=defects[0],parseval_relative_defect=defects[1]))
            checks.append(dict(**hashes,interface_jump=fingerprint(trace)))
        arrays[f'q{order}_element_error_squared_integrals']=np.array(radial)
        arrays[f'q{order}_cell_mode_error_squared_integrals']=np.array(spectra)
        arrays[f'q{order}_interface_jump_squared_integrals']=np.array(traces)
        quadrature.append(dict(order=order,metrics=metrics,physical_fingerprints=checks,
            original_density_fingerprints_reproduced=True))
    if any(sha(directory/p) != h for p,h in inputs.items()):
        raise ValueError('Original accepted prefix changed')
    return dict(directory=str(directory),spec=spec,prefix_steps=steps,
        input_raw_sha256=inputs,physical_source_sha256=source['source_sha256'],quadrature=quadrature),arrays


def run(density_directory,destination):
    density_directory=Path(density_directory);destination=Path(destination)
    if destination.exists():raise ValueError('Fresh localization destination required')
    source=json.loads((density_directory/'provenance.json').read_text())
    if any(sha(p)!=h for p,h in source['source_files'].items()):raise ValueError('Original density observer required')
    observed=json.loads((density_directory/'summary.json').read_text())
    if not observed['terminal']:raise ValueError('Complete independent density observation required')
    manifest_path=Path(observed['manifest'])
    if sha(manifest_path)!=observed['manifest_sha256']:raise ValueError('Original immutable prefix snapshots required')
    manifest=json.loads(manifest_path.read_text())
    if len(manifest['cases'])!=len(observed['cases']):raise ValueError('Complete matching cases required')
    inputs={str(density_directory/'summary.json'):sha(density_directory/'summary.json'),
        str(manifest_path):sha(manifest_path)}
    for case in observed['cases']:
        path=density_directory/case['arrays_file']
        if sha(path)!=case['arrays_sha256']:raise ValueError('Unchanged original density time series required')
        inputs[str(path)]=sha(path)
    identity=provenance(['examples.trajectory_density_localization'])
    destination.mkdir(parents=True);freeze(destination,identity);cases=[]
    for case,original in zip(manifest['cases'],observed['cases']):
        record,arrays=calculate(case,original);path=destination/f'case_{len(cases):02}.npz'
        np.savez_compressed(path,**arrays);record.update(arrays_file=path.name,arrays_sha256=sha(path));cases.append(record)
        print('density localization',record['directory'],flush=True)
    if any(sha(p)!=h for p,h in inputs.items()) or provenance(['examples.trajectory_density_localization'])!=identity:
        raise ValueError('Observation inputs or source changed')
    result=dict(terminal=True,source_sha256=identity['source_sha256'],density_directory=str(density_directory),
        input_sha256=inputs,cases=cases,
        scope='Full accepted-prefix density error partitions and actual internal trace jumps. Cell-translation modes fold aliases. Exact geometric partitions use field cells subdividing fluid elements. Surface jumps have distinct dimensions from volume errors and are normalized by background trace density. No cause, continuum stability, local conservation or full-duration accuracy is inferred.')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (destination/'resource_observation.json').write_text(json.dumps(dict(
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,free_bytes_after=shutil.disk_usage('.').free,
        scope='Own-process macOS cumulative peak and final disk observation; no continuous audit guard'),indent=2)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('density_directory');parser.add_argument('destination')
    args=parser.parse_args();run(args.density_directory,args.destination)
