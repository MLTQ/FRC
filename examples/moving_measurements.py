"""Recorded evolving free-material/vacuum verification with atomic partial records."""
from pathlib import Path
from time import perf_counter
import argparse
import hashlib
import json
import numpy as np
from .moving_cases import MovingCase, provenance, snapshot, errors
from .verification_records import write_record


def protocol():
    cases=[MovingCase(f"{kind}_{n}",kind,n,.1 if kind == "expansion" else .04)
           for kind in ("expansion","frc") for n in (4,8,16,32)]
    cases += [MovingCase(f"source_{n}","source",n,.04,.0005) for n in (4,8,16)]
    cases += [MovingCase(f"{kind}_time_{dt}",kind,8,.04,dt)
              for kind in ("expansion","frc","source") for dt in (.002,.001,.0005,.00025)]
    cases += [MovingCase("frc_far_matching","frc",16,.04,matching_radius=3.),
              MovingCase("frc_modes4","frc",16,.04,modes=4),
              MovingCase("frc_quadrature7","frc",16,.04,quadrature=7)]
    return cases


def measure(case, directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    source=provenance();start=perf_counter();s=case.build()
    row={"name":case.name,"completed":False,"configuration":case.configuration(),
         "provenance":source,"initial":s.diagnostics(),"initial_errors":errors(s,case.kind),
         "setup_seconds":perf_counter()-start,"history":[]}
    arrays={};saved=[];start=perf_counter();last=start
    partial=directory/f"{case.name}.progress.json"
    def progress(current,force=False):
        nonlocal last
        now=perf_counter()
        if force or now-last>=30:
            write_record(partial,{**row,"accepted_time_s":current.time,"accepted_steps":current.steps,
                                  "elapsed_seconds":now-start})
            print(case.name,"t=",current.time,"steps=",current.steps,flush=True);last=now
    progress(s,True)
    try:
        for i,end in enumerate(np.linspace(0.,case.end_time,11)):
            s.advance(float(end),on_step=progress)
            row["history"].append({**s.diagnostics(),"errors":errors(s,case.kind)})
            if i in (0,5,10):
                arrays.update({f"t{i}_{key}":value for key,value in snapshot(s).items()})
                saved.append({"prefix":f"t{i}","time_s":s.time})
        row["completed"]=True
    except (ValueError,RuntimeError) as error:
        row["failure"]=str(error)
    row.update(final=s.diagnostics(),final_errors=errors(s,case.kind),advance_seconds=perf_counter()-start)
    path=directory/f"{case.name}.npz";np.savez_compressed(path,**arrays)
    row["fields"]={"file":path.name,"sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"times":saved,
                   "units":{"coordinates":"m, components r/z","velocity":"m/s, components r/z",
                            "induced_a":"T, a=psi/r^2; vacuum modes last","total_material_a":"T",
                            "flux_labels":"initial a in T","node_mass":"kg","cell_mass":"kg",
                            "pressure_invariant":"Pa m^(3 gamma)","density":"kg/m^3","pressure":"Pa",
                            "volume":"m^3","bias":"T"}}
    if provenance()!=source:
        raise RuntimeError("Source changed during this measurement; do not label it reproducible")
    write_record(directory/f"{case.name}.json",row);partial.unlink(missing_ok=True)
    print("Finished",case.name,"completed=",row["completed"],flush=True)
    return row


def run(directory="outputs/moving_verification/cases",names=None,reuse=False):
    source=provenance()
    for case in protocol():
        if names and case.name not in names:
            continue
        path=Path(directory)/f"{case.name}.json"
        if reuse and path.exists():
            old=json.loads(path.read_text())
            if old["completed"] and old["provenance"]==source and old["configuration"]==case.configuration():
                print("Reuse",case.name,flush=True);continue
        measure(case,directory)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--names",nargs="*")
    parser.add_argument("--reuse",action="store_true");args=parser.parse_args()
    run(names=args.names,reuse=args.reuse)
