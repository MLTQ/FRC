"""Explicit experimental first-order boundary strip; production defaults unchanged."""
from dataclasses import asdict
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from plasma_bench.conservative_benchmarks import from_split_initial
from plasma_bench.shock_benchmarks import standing_shock
from .conservative_measurements import provenance,shock_errors
from .verification_records import write_record


def install_probe(solver):
    """Select baseline fluxes incident to wall/end cells before admissibility.

    The fixed strip is explicit spatial order reduction, distinct from later
    adaptive fallback. Its zero adaptive-limiter count does not mean full order 2.
    """
    if hasattr(solver,"boundary_order_probe"):
        raise ValueError("Boundary-order probe is already installed")
    if solver.parameters.spatial_order!=2:
        raise ValueError("Probe requires a second-order interior baseline")
    forced=np.zeros((solver.mesh.nr,solver.mesh.nz),dtype=bool)
    forced[-1,:]=True;forced[:,[0,-1]]=True
    original=solver._fluxes
    def fluxes(state,time,spatial_order=None):
        order=solver.parameters.spatial_order if spatial_order is None else spatial_order
        high,powers=original(state,time,order)
        if order==1:
            return high,powers
        low,_=original(state,time,1)
        mixed,_=high.local_fallback(low,forced)
        return mixed,powers
    affected=forced.copy()
    for di in (-1,0,1):
        for dj in (-1,0,1):
            i0,i1=max(0,-di),min(solver.mesh.nr,solver.mesh.nr-di)
            j0,j1=max(0,-dj),min(solver.mesh.nz,solver.mesh.nz-dj)
            affected[i0:i1,j0:j1] |= forced[i0+di:i1+di,j0+dj:j1+dj]
    solver._fluxes=fluxes
    solver.boundary_order_probe={"description":"First-order shared faces/corners and metric source incident to one wall/end cell strip; second-order interior",
                                 "forced_cell_count":int(forced.sum()),"affected_cell_count":int(affected.sum()),
                                 "adaptive_limiter_counts_exclude_this_fixed_order_reduction":True}
    return solver


def run(cells=96,cfl=.2,directory="outputs/compatible_energy_verification/boundary_order_probe"):
    split,reference=standing_shock(cells,cfl=cfl)
    solver=install_probe(from_split_initial(split,energy_quadrature="staggered",boundary_poynting="variational",time_order=3))
    source=provenance();source["boundary_probe_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    record={"experimental":True,"completed":False,"radial_cells":cells,"cfl":cfl,
            "parameters":asdict(solver.parameters),"fluid_mesh":asdict(solver.mesh),
            "field_grid":asdict(solver.mapping.field),"reference":asdict(reference),
            "bias_t":solver.sources.bias,"coils":[],"provenance":source,
            "boundary_order_probe":solver.boundary_order_probe,"history":[solver.diagnostics()]}
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    path=directory/f"shock_{cells}_cfl{cfl}.json"
    try:
        for end in (.005,.01,.015,.02):
            solver.advance(end);record["history"].append(solver.diagnostics())
            write_record(path,record)
            print(cells,cfl,end,record["history"][-1]["response_energy_balance_residual_j"],flush=True)
        record["completed"]=True
    except (ValueError,RuntimeError) as error:
        record["failure"]=str(error)
    record["diagnostics"]=solver.diagnostics()
    record["errors"],record["profile"]=shock_errors(solver,reference)
    write_record(path,record)
    return record


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells",type=int,default=96)
    parser.add_argument("--cfl",type=float,default=.2)
    args=parser.parse_args();run(args.cells,args.cfl)
