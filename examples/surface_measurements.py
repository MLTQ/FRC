"""Surface-normal, angular and temporal resolution at unchanged physical eta."""
from pathlib import Path
import argparse
import hashlib
import json
from .resistive_records import ResistiveCase,provenance as base_provenance
from .resistive_dynamic_measurements import measure
from .surface_references import run as references


def provenance():
    result=base_provenance()
    for name in ("surface_measurements","surface_references"):
        path=Path(f"examples/{name}.py")
        result["files"][str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    result["source_sha256"]=hashlib.sha256("\n".join(f"{p}:{v}" for p,v in result["files"].items()).encode()).hexdigest()
    return result


def protocol():
    cases=[ResistiveCase(f"source_g{g:g}_{n}","source",n,radial_grading=g)
           for g in (1.,1.5,2.) for n in (8,16,32)]
    cases += [ResistiveCase(f"frc_g{g:g}_{n}","frc",n,radial_grading=g)
              for g in (1.5,2.) for n in (16,32)]
    cases += [ResistiveCase("source_g2_16_angular","source",16,radial_grading=2.,angular_cells=64),
              ResistiveCase("source_g2_32_angular","source",32,radial_grading=2.,angular_cells=128),
              ResistiveCase("source_g2_16_time","source",16,radial_grading=2.,max_timestep=.00025),
              ResistiveCase("source_g2_32_cfl","source",32,radial_grading=2.,cfl=.075,max_timestep=.00025),
              ResistiveCase("source_g2_16_tolerance","source",16,radial_grading=2.,diffusion_tolerance=1e-5),
              ResistiveCase("source_g2_16_matching","source",16,radial_grading=2.,matching_radius=3.),
              ResistiveCase("source_g2_16_quadrature","source",16,radial_grading=2.,quadrature=7)]
    return cases


def run(directory="outputs/surface_verification",names=None,reuse=False,include_references=False):
    source=provenance();directory=Path(directory)
    if include_references:
        references(directory/"references",source,provenance)
        return
    known={c.name for c in protocol()}
    if names and not set(names)<=known:
        raise ValueError("Unknown surface study case")
    for case in protocol():
        if names and case.name not in names:
            continue
        path=directory/"dynamic"/f"{case.name}.json"
        if reuse and path.exists():
            old=json.loads(path.read_text());fields=path.parent/old["fields"]["file"]
            if (old["completed"] and old["provenance"]==source and old["configuration"]==case.configuration()
                    and fields.exists() and hashlib.sha256(fields.read_bytes()).hexdigest()==old["fields"]["sha256"]):
                print("Reuse",case.name,flush=True);continue
        measure(case,directory/"dynamic",fingerprint=provenance)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--names",nargs="*")
    parser.add_argument("--reuse",action="store_true");parser.add_argument("--references",action="store_true")
    args=parser.parse_args();run(names=args.names,reuse=args.reuse,include_references=args.references)
