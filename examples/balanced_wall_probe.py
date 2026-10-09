"""Replay an exact rejected state with fully traced dual-residual refinement."""
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.balanced_wall_plasma import BalancedWallPlasma
from plasma_bench.regular_plasma_step import field_from_packed
from plasma_bench.regular_current import RegularCurrent
from .archived_contact_mode import ArchivedContactMode
from .radial_contact_fixture import fixture
from .implicit_plasma_records import provenance, freeze, sha


class TracedPlasma(BalancedWallPlasma):
    def operators(self, model):
        result = super().operators(model)
        self.last_induction = result[-1]
        return result


def run(destination, original='outputs/radial16_fine_rates', index=1):
    destination=Path(destination); original=Path(original)
    if destination.exists(): raise ValueError('Refusing overwrite')
    metadata=json.loads((original/f'case_{index:02}.json').read_text())
    path=original/metadata['arrays_file']
    if sha(path)!=metadata['arrays_sha256']:raise ValueError('Captured input hash mismatch')
    saved=dict(np.load(path)); mode=ArchivedContactMode(); spec=metadata['spec']
    identity=provenance(['examples.balanced_wall_probe']); destination.mkdir(parents=True);freeze(destination,identity)
    start=time.monotonic();d,u,s,bias,du,db=fixture(spec,mode)
    if not all(np.array_equal(saved[k],v) for k,v in [('state',u),('fluid_perturbation',du),('magnetic_perturbation',db)]):
        raise ValueError('Regenerated fixture differs from captured input')
    s=TracedPlasma(s.mapping,s.ion_mass,cache_geometry=s.cache_geometry,**s.solver_options)
    state=saved['minus_state'].copy(); packed=saved['minus_magnetic'].copy()
    field=field_from_packed(s.geometry,packed,bias.axial_bias)
    arrays=dict(state=state.copy(),magnetic=packed.copy())
    record=dict(source_sha256=identity['source_sha256'],input_sha256=sha(path),input_directory=str(original),spec=spec)
    try:
        rate=s.derivative(state,field)
        arrays.update(fluid_rate=rate.fluid,magnetic_rate=rate.magnetic,current_rate=rate.current,electric=rate.electric)
        model=RegularCurrent(s.mapping,state[...,0],s.ion_mass);j=model.current(packed)
        powers=[np.sum(d.volume[...,None]*d.fields(state)['energy_gradient']*rate.fluid),
                j@(model.inertia@rate.current)+np.sum(model.density_energy_gradient(j)*rate.fluid[...,0]),
                packed@(s.geometry.hessian@rate.magnetic)]
        record.update(status='complete',residual=rate.residual,diagnostics=rate.diagnostics,
            independent_powers=[float(v) for v in powers],relative_power_defect=float(abs(sum(powers))/sum(abs(v) for v in powers)))
    except (RuntimeError,ValueError) as error:
        record.update(status='failed',failure_type=type(error).__name__,failure_message=str(error))
    record['trace']=getattr(s.last_induction,'attempt_trace',[])
    record['unchanged_input']=bool(np.array_equal(state,saved['minus_state']) and np.array_equal(field.packed_increment(s.geometry),saved['minus_magnetic']))
    if not record['unchanged_input']:raise RuntimeError('Input mutation')
    if provenance(['examples.balanced_wall_probe'])!=identity:raise RuntimeError('Source changed')
    path=destination/'probe.npz';np.savez_compressed(path,**arrays)
    record.update(terminal=True,arrays_file=path.name,arrays_sha256=sha(path),elapsed_seconds=time.monotonic()-start)
    (destination/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    print(record['status'],record.get('residual'),record.get('relative_power_defect'),record['trace'],flush=True)
    return record
