"""Cache costly overlap integrations using exact input and analysis fingerprints."""
from pathlib import Path
import argparse
import hashlib
import json
from .surface_analysis import compare,comparison_pairs,read_study
from .verification_records import write_record


def analysis_sources():
    return {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(f"examples/{n}.py") for n in ("surface_comparisons","surface_analysis","surface_sampling","material_overlap")]}


def measurements(directory,rows,arrays):
    directory=Path(directory)/"comparisons";directory.mkdir(parents=True,exist_ok=True)
    sources=analysis_sources();results=[]
    for test,reference in comparison_pairs():
        if test not in rows or reference not in rows:continue
        identity={"analysis_sources":sources,"test":rows[test],"reference":rows[reference],"orders":[4,6]}
        name=f"{test}__{reference}";path=directory/f"{name}.json"
        if path.exists():
            old=json.loads(path.read_text())
            if old.get("completed") and old["identity"]==identity:
                results.append(old);print("Reuse comparison",name,flush=True);continue
        low=compare(rows[test],rows[reference],arrays,4)
        high=compare(rows[test],rows[reference],arrays,6)
        if analysis_sources()!=sources:raise RuntimeError("Comparison source changed during integration")
        result={"completed":True,"identity":identity,"order4":low,"order6":high}
        write_record(path,result);results.append(result);print("Finished comparison",name,flush=True)
    return results


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--base",action="store_true");args=parser.parse_args()
    directory="outputs/surface_verification"
    _,rows,arrays=read_study(directory,include_extra=not args.base)
    measurements(directory,rows,arrays)
