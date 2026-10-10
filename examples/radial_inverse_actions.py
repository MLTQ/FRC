"""Decompose the frozen opposing residual gap into mass-inverse actions."""
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.wall_induction_core import WallInductionCore
from .archived_contact_mode import ArchivedContactMode
from .radial_contact_fixture import fixture
from .implicit_plasma_records import provenance,freeze,sha


def run(destination,original='outputs/radial16_roundoff_probe'):
    destination=Path(destination);root=Path(original);record=json.loads((root/'summary.json').read_text());path=root/record['arrays_file']
    if destination.exists():raise ValueError('Refusing overwrite')
    if not record['terminal'] or sha(path)!=record['arrays_sha256']:raise ValueError('Exact terminal probe required')
    raw=dict(np.load(path));identity=provenance(['examples.radial_inverse_actions']);destination.mkdir(parents=True);freeze(destination,identity);start=time.monotonic()
    _,_,s,_,_,_=fixture(record['spec'],ArchivedContactMode());model=RegularCurrent(s.mapping,raw['state'][...,0],s.ion_mass)
    core=WallInductionCore(model,tolerance=record['tolerance']);geo=model.geometry;inverse=core.mass_inverse;f=raw['weak_remainder'];base=inverse(f)
    norm=np.linalg.norm;relative=lambda a,b:float(norm(a)/max(norm(b),np.finfo(float).tiny));records=[];arrays={}
    for i in range(len(record['trace'])):
        y=raw[f'refinement_{i}_solution'];rhs=raw[f'refinement_{i}_rhs'];rate=core.lift@y
        pre=core.dual_lift@y;post=geo.hessian@rate
        jp=inverse(geo.faraday.T@pre);jq=inverse(geo.faraday.T@post)
        wp=model.inertia@jp;wq=model.inertia@jq;cp=inverse(wp);cq=inverse(wq);combined=inverse(f-wq)
        split=base-cq;distributive=combined-base+cq
        mass_gap=core.dual_lift.T@(geo.faraday@distributive)
        current_gap=core.dual_lift.T@(geo.faraday@(cp-cq))
        metric_gap=core.lift.T@pre-core.dual_lift.T@rate
        observed=raw[f'refinement_{i}_operator_residual']+raw[f'refinement_{i}_physical_dual_residual']
        predicted=mass_gap+current_gap+metric_gap
        delta=geo.faraday@split-rate;physical=float(np.sqrt(delta@(geo.hessian@delta))/np.sqrt(rate@(geo.hessian@rate)))
        row=dict(refinement=i,mass_inverse_distributivity_relative_electric=relative(distributive,combined),
            observed_gap_over_rhs=relative(observed,rhs),mass_gap_over_rhs=relative(mass_gap,rhs),current_gap_over_rhs=relative(current_gap,rhs),
            metric_gap_over_rhs=relative(metric_gap,rhs),unexplained_fraction=relative(predicted-observed,observed),
            original_operator_residual=record['trace'][i]['operator'],split_faraday_dual=relative(core.dual_lift.T@delta,rhs),split_faraday_physical=physical,
            split_original_ohm_backward=relative(model.mass@split+wq-f,abs(f)+abs(wq)))
        records.append(row)
        for key,v in dict(observed_gap=observed,mass_gap=mass_gap,current_gap=current_gap,metric_gap=metric_gap,electric_distributivity=distributive,split_electric=split).items():arrays[f'refinement_{i}_{key}']=v
    path=destination/'actions.npz';np.savez_compressed(path,**arrays)
    if provenance(['examples.radial_inverse_actions'])!=identity:raise ValueError('Source changed')
    result=dict(terminal=True,source_sha256=identity['source_sha256'],input_arrays_sha256=record['arrays_sha256'],records=records,
        full_plasma_rate_accepted=False,arrays_file=path.name,arrays_sha256=sha(path),elapsed_seconds=time.monotonic()-start)
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print('inverse action diagnosis',records,flush=True)
