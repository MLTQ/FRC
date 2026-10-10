"""Original full physical term coefficients at archived signed states."""
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
    root=Path(directory);record=json.loads((root/'summary.json').read_text())['cases'][index]
    if record['status']!='complete':raise ValueError('Accepted archived case required')
    path=root/record['arrays_file']
    if sha(path)!=record['arrays_sha256']:raise ValueError('Original arrays changed')
    saved=dict(np.load(path));mode=ArchivedContactMode();d,u,s,bias,du,db=fixture(record['spec'],mode)
    if any(not np.array_equal(saved[k],v) for k,v in [('state',u),('fluid_perturbation',du),('magnetic_perturbation',db)]):raise ValueError('Different physical fixture')
    arrays={k:v for k,v in saved.items() if not k.startswith('q')};checks=[];operators=None
    for sign in ('minus','plus'):
        state=saved[sign+'_state'];packed=saved[sign+'_magnetic'];field=field_from_packed(s.geometry,packed,bias.axial_bias)
        model=RegularCurrent(s.mapping,state[...,0],s.ion_mass);geo=model.geometry
        bulk,_=s.fluid.derivative(state);em=model.electromagnetic_pair(state,field,bulk[...,0])
        if operators is None:operators=(CurrentTransport(model),RegularSpecies(model))
        transport,species=(copy(op) for op in operators)
        transport.model=model;species.model=model
        tr=transport.pair(state,em['current'],bulk[...,0]);sp=species.pair(state,em['current'])
        fluid=dict(bulk=bulk,lorentz=em['fluid_rate'],drift=tr['fluid_rate'],species=sp['fluid_rate'])
        total=bulk+em['fluid_rate']+tr['fluid_rate']+sp['fluid_rate'];reaction=np.zeros_like(total)
        reaction[-1,-1,:,1]=-total[-1,-1,:,1];fluid['wall_reaction']=reaction;total=total+reaction
        if not np.array_equal(total,saved[sign+'_fluid_rate']):raise ValueError('Full fluid RHS differs from archived accepted rate')
        inverse=lambda value:np.concatenate([solve(v) for solve,v in zip(geo.solves,geo.split_electric(value))])
        weak=dict(ideal=em['ideal_forcing'],hall=em['hall_forcing'],density=em['density_forcing'],transport=tr['forcing'],
            pressure_full=sp['forcing'],pressure_reduced=sp['reduced_forcing'],inertia=model.inertia@saved[sign+'_current_rate'])
        electric={key:inverse(value) for key,value in weak.items()};electric['pressure_potential']=sp['potential_electric']
        forcing=em['forcing']+tr['forcing']+sp['forcing'];actual=saved[sign+'_electric'];inertia=weak['inertia']
        scale=np.linalg.norm(abs(model.mass@actual)+abs(inertia)+abs(forcing))
        ohm=float(np.linalg.norm(model.mass@actual-inertia-forcing)/scale)
        if ohm>s.solver_options['tolerance']:raise ValueError('Archived full Ohm residual rejected')
        reconstructed=sum(electric[k] for k in ('ideal','hall','density','transport','pressure_full','inertia'))
        checks.append(dict(sign=sign,fluid_bitwise_equal=True,original_full_ohm_backward=ohm,
            electric_sum_relative=float(np.linalg.norm(reconstructed-actual)/np.linalg.norm(actual))))
        for group,values in [('fluid',fluid),('electric',electric),('weak',weak)]:
            arrays.update({f'{sign}_{group}_{k}':v for k,v in values.items()})
    return dict(input_directory=str(root),input_index=index,input_arrays_sha256=sha(path),spec=record['spec'],epsilon=record['epsilon'],
        reference_sha256=mode.input_sha256,checks=checks),arrays


def run(destination,cases):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    identity=provenance(['examples.radial_term_capture']);destination.mkdir(parents=True);freeze(destination,identity);records=[]
    for directory,index in cases:
        start=time.monotonic();record,arrays=calculate(directory,index);path=destination/f'case_{len(records):02}.npz'
        np.savez_compressed(path,**arrays);record.update(arrays_file=path.name,arrays_sha256=sha(path),elapsed_seconds=time.monotonic()-start);records.append(record)
        print('captured terms',record['spec'],record['checks'],flush=True)
    if provenance(['examples.radial_term_capture'])!=identity:raise ValueError('Source changed')
    (destination/'summary.json').write_text(json.dumps(dict(terminal=True,source_sha256=identity['source_sha256'],cases=records),indent=2)+'\n')
