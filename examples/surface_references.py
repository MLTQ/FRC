"""Independent gas expansion and magnetic decay on interface-refined meshes."""
from pathlib import Path
import hashlib
import numpy as np
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.moving_plasma import MovingPlasma
from .moving_cases import errors,snapshot
from .resistive_decay_measurements import decay_case
from .verification_records import write_record


def expansion_case(n,grading,end_time=.04):
    d=MovingDomain(n,2*n,radial_grading=grading)
    s=MovingPlasma(d,0.,pressure=lambda r,z:1-r*r-z*z,max_timestep=.0005)
    row={"name":f"expansion_g{grading:g}_{n}","completed":False,
         "configuration":{"radial_cells":n,"angular_cells":2*n,"radial_grading":grading,
          "radius_m":1.,"matching_radius_m":2.,"vacuum_layers":n,"modes":8,"quadrature":5,
          "gamma":5/3,"density_kg_m3":1.,"pressure_pa":"1-r^2-z^2","initial_velocity_m_s":0.,
          "physical_resistivity_ohm_m":0.,"initial_a_t":0.,"end_time_s":end_time,
          "max_timestep_s":.0005,"cfl":.15,"reference":"lambda=sqrt(1+2t^2); no heat source"},
         "initial":s.diagnostics()}
    arrays={f"initial_{k}":v for k,v in snapshot(s).items()}
    s.advance(end_time)
    row.update(completed=True,final=s.diagnostics(),errors=errors(s,"expansion"))
    arrays.update({f"final_{k}":v for k,v in snapshot(s).items()})
    return row,arrays


def run(directory,source,fingerprint):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);rows=[]
    cases=[(f"decay_g{g:g}_{n}",lambda n=n,g=g:decay_case(f"decay_g{g:g}_{n}",n,radial_grading=g))
           for g in (1.5,2.) for n in (8,16,32)]
    cases += [(f"expansion_g{g:g}_{n}",lambda n=n,g=g:expansion_case(n,g))
              for g in (1.5,2.) for n in (4,8,16)]
    cases += [("decay_g2_16_angular",lambda:decay_case("decay_g2_16_angular",16,radial_grading=2.,angular_cells=64))]
    for name,compute in cases:
        row,arrays=compute();path=directory/f"{name}.npz";np.savez_compressed(path,**arrays)
        row.update(provenance=source,fields={"file":path.name,"sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
                   "units":"SI; decay component matches resistive_decay_measurements; gas component matches moving_cases"})
        if fingerprint()!=source:
            raise RuntimeError("Source changed during surface-reference measurements")
        write_record(directory/f"{name}.json",row);rows.append(row)
        print("Finished surface reference",name,flush=True)
    write_record(directory/"metrics.json",{"completed":True,"provenance":source,"cases":rows})
