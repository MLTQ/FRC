"""Explicit recovery of missing case files after a host storage failure."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
from pathlib import Path
from .bulk_thermal_evolution import protocol,run_case
from .bulk_thermal_records import provenance,sha


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    parser.add_argument('--indices',type=int,nargs='*',default=[]);parser.add_argument('--assemble',action='store_true')
    args=parser.parse_args();directory=args.directory;identity=provenance()
    if json.loads((directory/'provenance.json').read_text())['source_sha256']!=identity['source_sha256']:
        raise RuntimeError('Cannot recover with a changed scientific source')
    specs=protocol()
    if len(set(args.indices))!=len(args.indices) or any(i<0 or i>=len(specs) for i in args.indices):raise ValueError('Distinct protocol indices required')
    for index in args.indices:
        if list(directory.glob(f'transport_{index:02d}*')):raise RuntimeError('Preserve existing evidence before requesting a recomputation')
    if args.indices:
        note=directory/'disk_recovery.json'
        if note.exists():raise RuntimeError('Do not overwrite a previous recovery record')
        record=dict(source_sha256=identity['source_sha256'],indices=args.indices,helper_source_sha256=sha(Path(__file__)),
                    reason='Host disk filled during the original run; only explicitly missing files are recomputed, no numerical step retry or state repair',cases=[])
        note.write_text(json.dumps(record,indent=2)+'\n')
        with ProcessPoolExecutor(max_workers=2) as pool:
            futures=[pool.submit(run_case,i,specs[i],directory,identity['source_sha256']) for i in args.indices]
            for future in as_completed(futures):
                index,result=future.result();record['cases'].append(dict(index=index,status=result['status']))
                note.write_text(json.dumps(record,indent=2)+'\n')
    if args.assemble:
        records=[]
        for index,spec in enumerate(specs):
            r=json.loads((directory/f'transport_{index:02d}.json').read_text())
            if r['spec']!=spec or r['status']!='complete' or r['source_sha256']!=identity['source_sha256'] or sha(directory/r['arrays_file'])!=r['arrays_sha256']:
                raise RuntimeError('Complete unchanged cases required for assembly')
            records.append(r)
        output=directory/'transports.json'
        if output.exists():raise RuntimeError('Do not overwrite an existing protocol')
        output.write_text(json.dumps(records,indent=2)+'\n')
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Source changed during recovery')


if __name__=='__main__':main()
