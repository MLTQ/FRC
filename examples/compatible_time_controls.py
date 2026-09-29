"""First-order spatial controls isolate fine-boundary temporal sensitivity."""
from pathlib import Path
import hashlib
import json
from .conservative_measurements import measure,provenance
from .verification_records import write_record


def run(directory="outputs/compatible_energy_verification/temporal_controls",*,reuse=False):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    source=provenance()
    source["temporal_controls_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    rows=[]
    for cfl in (.2,.1):
        name=f"first_order_192_cfl{cfl}"
        path=directory/f"{name}.json"
        options={"cfl":cfl,"spatial_order":1,"energy_quadrature":"staggered",
                 "boundary_poynting":"variational","time_order":3}
        cached=json.loads(path.read_text()) if reuse and path.exists() else None
        if cached and cached.get("provenance")==source and cached.get("input_options")==options:
            row=cached
        else:
            def progress(partial):
                partial.update(name=name,provenance=source,input_options=options)
                write_record(directory/f"{name}.progress.json",partial)
                d=partial["diagnostics"]
                print(name,d["time_s"],d["steps"],d["response_energy_balance_residual_j"],flush=True)
            row=measure("shock",192,progress=progress,**options)
            row.update(name=name,provenance=source,input_options=options)
            write_record(path,row)
            (directory/f"{name}.progress.json").unlink(missing_ok=True)
        rows.append(row)
    return rows


if __name__=="__main__":
    run(reuse=True)
