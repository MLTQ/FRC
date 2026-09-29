"""Preserve uniform-bias equilibrium precision controls without hiding rejection."""
from pathlib import Path
import hashlib
import numpy as np
from scipy.constants import atomic_mass,elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.current_thermal import CurrentThermal
from plasma_bench.current_thermal_midpoint import CurrentThermalMidpoint
from plasma_bench.hall_midpoint import HallConvergenceError
from .current_thermal_records import provenance
from .verification_records import write_record


def run(directory=Path('outputs/current_bias_failure')):
    directory.mkdir(parents=True,exist_ok=True);source=provenance();records=[]
    for nr in (8,16):
        for dt in (1e-8,1e-6):
            mesh=MHDMesh(.035,.12,nr,2*nr);n=np.full((nr+1,2*nr),1e20)
            model=CurrentThermal(mesh,n,39.948*atomic_mass)
            psi=.5*.005*mesh.corners[0]**2;g=np.zeros((nr,2*nr))
            ue=model.electron.from_pressure(n*5*E);ui=model.ion.from_pressure(n*E)
            initial=[v.copy() for v in (psi,g,ue,ui)]
            current=model.geometry.sample_electric(model.current(model.geometry.pack(psi,g)))
            record=dict(nr=nr,nz=2*nr,dt_s=dt,radius_m=.035,length_m=.24,bias_t=.005,density_m3=1e20,
                electron_ev=5.,ion_ev=1.,ion_mass=39.948*atomic_mass,tolerance=1e-10,max_iterations=80,
                exact_continuum_current_a_m2=0.,initial_weak_current_max_a_m2=float(np.linalg.norm(current,axis=-1).max()))
            arrays=dict(initial_psi=psi,initial_g=g,initial_electron_entropy=ue,initial_ion_entropy=ui)
            try:
                result=CurrentThermalMidpoint(model).step(psi,g,ue,ui,dt)
                record.update(completed=True,accepted_time_s=dt,nonlinear_residual=result.nonlinear_residual)
                arrays.update(final_psi=result.psi,final_g=result.g,final_electron_entropy=result.electron_entropy,final_ion_entropy=result.ion_entropy)
            except HallConvergenceError as error:
                record.update(completed=False,accepted_time_s=0.,reason=str(error),input_unchanged=all(np.array_equal(a,b) for a,b in zip(initial,(psi,g,ue,ui))))
                assert record['input_unchanged']
            path=directory/f'bias_{nr}_{dt:g}.npz';np.savez_compressed(path,**arrays)
            record['fields']=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest());records.append(record)
            print(record,flush=True)
    assert source==provenance()
    write_record(directory/'measurements.json',dict(observation_finished=True,provenance=source,cases=records,
        helper_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limitation='Uniform physical equilibrium exposes weak-current cancellation and returned-state inertia residual sensitivity; no accepted physical evolution is claimed for rejections.'))


if __name__=='__main__':run()
