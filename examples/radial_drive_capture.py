"""Physical driving terms on valid inputs, including rejected induction cases."""
from copy import copy
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.current_transport import CurrentTransport
from plasma_bench.regular_species import RegularSpecies
from plasma_bench.regular_plasma_step import field_from_packed
from .archived_contact_mode import ArchivedContactMode
from .radial_contact_fixture import fixture
from .implicit_plasma_records import provenance,freeze,sha


def calculate(directory,index):
    root=Path(directory);record=json.loads((root/f'case_{index:02}.json').read_text());path=root/record['arrays_file']
    if sha(path)!=record['arrays_sha256']:raise ValueError('Input archive changed')
    saved=dict(np.load(path));mode=ArchivedContactMode();d,u,s,bias,du,db=fixture(record['spec'],mode);epsilon=record['epsilon']
    arrays=dict(state=u,fluid_perturbation=du,magnetic_perturbation=db)
    if any(not np.array_equal(saved[k],v) for k,v in arrays.items()):raise ValueError('Fixture mismatch')
    checks=[];operators=None
    for sign,sgn in [('minus',-1),('plus',1)]:
        state=u+sgn*epsilon*du;packed=sgn*epsilon*db;field=field_from_packed(s.geometry,packed,bias.axial_bias)
        if sign+'_state' in saved and not np.array_equal(state,saved[sign+'_state']):raise ValueError('Signed state mismatch')
        arrays.update({sign+'_state':state,sign+'_magnetic':field.packed_increment(s.geometry)})
        model=RegularCurrent(s.mapping,state[...,0],s.ion_mass);geo=model.geometry;bulk,_=s.fluid.derivative(state)
        em=model.electromagnetic_pair(state,field,bulk[...,0])
        if operators is None:operators=(CurrentTransport(model),RegularSpecies(model))
        tr_op,sp_op=(copy(op) for op in operators);tr_op.model=model;sp_op.model=model
        tr=tr_op.pair(state,em['current'],bulk[...,0]);sp=sp_op.pair(state,em['current'])
        fluid=dict(bulk=bulk,lorentz=em['fluid_rate'],drift=tr['fluid_rate'],species=sp['fluid_rate'])
        total=bulk+em['fluid_rate']+tr['fluid_rate']+sp['fluid_rate'];reaction=np.zeros_like(total)
        reaction[-1,-1,:,1]=-total[-1,-1,:,1];fluid['wall_reaction']=reaction;fluid['total']=total+reaction
        equal=None
        if sign+'_fluid_rate' in saved:
            equal=bool(np.array_equal(fluid['total'],saved[sign+'_fluid_rate']))
            if not equal:raise ValueError('Saved full fluid rate mismatch')
        weak=dict(ideal=em['ideal_forcing'],hall=em['hall_forcing'],density=em['density_forcing'],transport=tr['forcing'],pressure_full=sp['forcing'],pressure_reduced=sp['reduced_forcing'])
        inverse=lambda v:np.concatenate([solve(part) for solve,part in zip(geo.solves,geo.split_electric(v))])
        electric={k:inverse(v) for k,v in weak.items()};electric['pressure_potential']=sp['potential_electric']
        for group,data in [('fluid',fluid),('weak',weak),('electric',electric)]:arrays.update({f'{sign}_{group}_{k}':v for k,v in data.items()})
        checks.append(dict(sign=sign,archived_full_fluid_bitwise_equal=equal))
    return dict(input_directory=str(root),input_index=index,input_arrays_sha256=sha(path),spec=record['spec'],epsilon=epsilon,
        original_full_plasma_status=record['status'],full_plasma_rate_accepted=False,reference_sha256=mode.input_sha256,checks=checks),arrays


def run(destination,cases):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    identity=provenance(['examples.radial_drive_capture']);destination.mkdir(parents=True);freeze(destination,identity);records=[]
    for directory,index in cases:
        start=time.monotonic();record,arrays=calculate(directory,index);path=destination/f'case_{len(records):02}.npz';np.savez_compressed(path,**arrays)
        record.update(arrays_file=path.name,arrays_sha256=sha(path),elapsed_seconds=time.monotonic()-start);records.append(record)
        (destination/f'case_{len(records)-1:02}.json').write_text(json.dumps(record,indent=2)+'\n')
        print('drive capture',record['spec'],record['original_full_plasma_status'],record['checks'],flush=True)
    if provenance(['examples.radial_drive_capture'])!=identity:raise ValueError('Source changed')
    (destination/'summary.json').write_text(json.dumps(dict(terminal=True,source_sha256=identity['source_sha256'],cases=records),indent=2)+'\n')
