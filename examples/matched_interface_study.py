"""Frozen original cell/interface budgets on independently refined physical modes."""
import json
from pathlib import Path
import time
import numpy as np
from .matched_mode_fixture import fixture
from .physical_cell_budget import measure
from .implicit_plasma_records import provenance,freeze,sha
from plasma_bench.regular_plasma_step import field_from_packed
from plasma_bench.regular_current import RegularCurrent


def calculate(spec):
    start=time.monotonic();mode,d,base,system,bias,du,db=fixture(spec)
    epsilon=spec.get('epsilon',1e-4);state=base+epsilon*du
    field=field_from_packed(system.geometry,epsilon*db,bias.axial_bias)
    rate=system.derivative(state,field)
    current=RegularCurrent(system.mapping,state[...,0],system.ion_mass).current(epsilon*db)
    arrays=dict(state=state,magnetic=epsilon*db,current=current,fluid_rate=rate.fluid,
        magnetic_rate=rate.magnetic,current_rate=rate.current,electric=rate.electric)
    quadrature=[]
    for order in (16,20):
        record,raw=measure(d,system,state,field,rate,order);quadrature.append(record)
        arrays.update({f'q{order}_{name}':value for name,value in raw.items()})
    return dict(spec=spec,epsilon=epsilon,input_sha256=mode.input_sha256,quadrature=quadrature,
        induction_residual=rate.residual,full_rhs_diagnostics=rate.diagnostics,
        elapsed_seconds=time.monotonic()-start),arrays


def run(destination,specs):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    seeds=['examples.matched_interface_study'];identity=provenance(seeds)
    destination.mkdir(parents=True);freeze(destination,identity);cases=[]
    for spec in specs:
        record,arrays=calculate(spec);raw=destination/f'case_{len(cases):02}.npz'
        np.savez_compressed(raw,**arrays);record.update(arrays_file=raw.name,arrays_sha256=sha(raw));cases.append(record)
        (destination/f'case_{len(cases)-1:02}.json').write_text(json.dumps(record,indent=2)+'\n')
        print('interface',spec,{k:v['cell_balance_relative_L2'] for k,v in record['quadrature'][1]['channels'].items()},flush=True)
    if provenance(seeds)!=identity:raise RuntimeError('Source changed')
    result=dict(terminal=True,source_sha256=identity['source_sha256'],cases=cases)
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
