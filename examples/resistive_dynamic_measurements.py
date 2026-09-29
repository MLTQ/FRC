"""Actual moving resistive plasma: mesh, source, time and boundary controls."""
from pathlib import Path
from time import perf_counter
import argparse
import hashlib
import json
import numpy as np
from .resistive_records import ResistiveCase,provenance,snapshot
from .verification_records import write_record


def protocol():
    cases=[ResistiveCase(f"{kind}_{n}",kind,n) for kind in ("frc","source") for n in (4,8,16,32)]
    cases += [ResistiveCase(f"source_wide_{n}","source",n,diffusivity=1.) for n in (16,32)]
    cases += [ResistiveCase(f"{kind}_time_{dt}",kind,8,max_timestep=dt)
              for kind in ("frc","source") for dt in (.001,.00025)]
    cases += [ResistiveCase(f"{kind}_far_matching",kind,16,matching_radius=3.) for kind in ("frc","source")]
    cases += [ResistiveCase(f"{kind}_tolerance",kind,16,diffusion_tolerance=1e-5) for kind in ("frc","source")]
    cases += [ResistiveCase("source_modes4","source",16,modes=4),
              ResistiveCase("source_quadrature7","source",16,quadrature=7)]
    return cases


def measure(case,directory,*,fingerprint=provenance):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    source=fingerprint();start=perf_counter();s=case.build()
    row={"name":case.name,"completed":False,"provenance":source,"configuration":case.configuration(),
         "initial":s.diagnostics(),"setup_seconds":perf_counter()-start,"history":[]}
    start=perf_counter();last=start;arrays={};times=[]
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
            row["history"].append(s.diagnostics())
            if i in (0,5,10):
                arrays.update({f"t{i}_{k}":v for k,v in snapshot(s).items()})
                times.append({"prefix":f"t{i}","time_s":s.time})
        row["completed"]=True
    except (ValueError,RuntimeError) as error:
        row["failure"]=str(error)
    row.update(final=s.diagnostics(),advance_seconds=perf_counter()-start)
    path=directory/f"{case.name}.npz";np.savez_compressed(path,**arrays)
    row["fields"]={"file":path.name,"sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"times":times,
        "units":{"coordinates":"m, r/z","velocity":"m/s, r/z","induced_a":"T, vacuum modes last",
                 "total_material_a":"T","flux_labels":"reference a in T, physically nonideal labels evolve",
                 "pressure":"Pa","density":"kg/m^3","volume":"m^3","node_mass":"kg","cell_mass":"kg",
                 "pressure_invariant":"Pa m^(3 gamma)","current_over_radius":"A/m^3",
                 "current_quadrature":"A/m^2","quadrature_coordinates":"m, r/z","cell_ohmic_heat":"J","bias":"T"}}
    if fingerprint()!=source:
        raise RuntimeError("Source changed during resistive measurement")
    write_record(directory/f"{case.name}.json",row);partial.unlink(missing_ok=True)
    print("Finished",case.name,"completed=",row["completed"],flush=True)
    return row


def run(directory="outputs/resistive_moving_verification/dynamic",names=None,reuse=False):
    source=provenance()
    for case in protocol():
        if names and case.name not in names:
            continue
        path=Path(directory)/f"{case.name}.json"
        if reuse and path.exists():
            old=json.loads(path.read_text())
            if old["completed"] and old["configuration"]==case.configuration() and old["provenance"]==source:
                print("Reuse",case.name,flush=True);continue
        measure(case,directory)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--names",nargs="*")
    parser.add_argument("--reuse",action="store_true");args=parser.parse_args()
    run(names=args.names,reuse=args.reuse)
