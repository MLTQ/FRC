"""Run and strictly assemble the finite-mass current/thermal protocol."""
import argparse
from pathlib import Path
from time import perf_counter
import hashlib
import json
from .current_thermal_cases import protocol
from .current_thermal_evolution import measure
from .current_symbol_measurement import measure as symbol_measure
from .current_thermal_records import provenance,save
from .verification_records import write_record


def assemble(directory):
    source=provenance();records=[]
    for case in protocol():
        record=json.loads((directory/(case['name']+'.json')).read_text())
        assert record['configuration']==case and record['provenance']==source
        assert record['status'] in ('succeeded','rejected') and record['completed']==(record['status']=='succeeded')
        assert hashlib.sha256((directory/record['fields']['file']).read_bytes()).hexdigest()==record['fields']['sha256']
        records.append(record)
    symbol=json.loads((directory/'principal_symbol.json').read_text());assert symbol['provenance']==source
    write_record(directory/'measurements.json',dict(protocol_finished=True,provenance=source,cases=records,principal_symbol=symbol))
    print('Assembled',len(records),'records;',sum(r['completed'] for r in records),'completed',source['source_sha256'],flush=True)


def run(directory,group=None):
    if group not in (None,'original','wave','symbol') and group not in {c['name'] for c in protocol()}:
        raise ValueError('Unknown current/thermal measurement group or case')
    directory.mkdir(parents=True,exist_ok=True);source=provenance()
    if group in (None,'symbol'):
        write_record(directory/'principal_symbol.json',dict(provenance=source,**symbol_measure()))
    for case in protocol():
        if group is not None and case['kind']!=group and case['name']!=group:continue
        start=perf_counter();record,arrays=measure(case)
        save(directory,case,source,record|dict(elapsed_seconds=perf_counter()-start),arrays)
        print(case['name'],record['status'],record['metrics'],record['failure'],flush=True)
    if group is None:assemble(directory)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',default='outputs/current_thermal_verification')
    parser.add_argument('--group');parser.add_argument('--assemble',action='store_true');args=parser.parse_args()
    if args.assemble:assemble(Path(args.directory))
    else:run(Path(args.directory),args.group)
