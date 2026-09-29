"""Frozen rate protocol for complete moving conservative finite-mass dynamics."""
from dataclasses import asdict,replace
from pathlib import Path
import argparse
import hashlib
from time import perf_counter
import numpy as np
from scipy.constants import elementary_charge as E,atomic_mass
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.finite_mass_rhs import rhs
from .moving_species_reference import MovingSpeciesReference
from .finite_mass_wave_reference import FiniteMassWave
from .finite_mass_manufactured import initialize,ManufacturedSources
from .finite_mass_wave_case import initialize as wave_initial
from .finite_mass_rate_norms import measure
from .finite_mass_records import provenance
from .verification_records import write_record


def protocol():
    ref=MovingSpeciesReference();cases=[]
    def add(name,nr,family,reference,phase=0.):
        cases.append(dict(name=name,nr=nr,family=family,reference=asdict(reference),phase=phase))
    for nr in (8,16,32,64,128):add(f'moving_{nr}',nr,'manufactured',ref,.23)
    for nr in (8,16,32,64):add(f'target_{nr}',nr,'manufactured',replace(ref,density=1e20),.23)
    base=FiniteMassWave(.035,.24,1e18,39.948*atomic_mass,.005,6e18*E,gamma_i=1.4,amplitude=1.)
    for branch in (0,1,2):
        for nr in (8,16,32):add(f'wave_{branch}_{nr}',nr,'wave',replace(base,branch=branch))
    add('wave_0_64',64,'wave',replace(base,branch=0))
    for branch in (0,2):add(f'weak_{branch}_32',32,'wave',replace(base,branch=branch,amplitude=.1))
    return cases


def case_state(case):
    ref=(MovingSpeciesReference if case['family']=='manufactured' else FiniteMassWave)(**case['reference'])
    mesh=MHDMesh(ref.radius,ref.length/2,case['nr'],2*case['nr'])
    if case['family']=='manufactured':
        time=case['phase']/ref.frequency;state,parameters=initialize(mesh,ref,time);source=ManufacturedSources(mesh,ref)
    else:time=0.;state,parameters=wave_initial(mesh,ref);source=None
    return ref,mesh,state,parameters,source,time


def run(directory,group=None):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);frozen=provenance()
    if group is not None and group not in {c['name'] for c in protocol()}:raise ValueError('Unknown case')
    for case in protocol():
        if group is not None and case['name']!=group:continue
        start=perf_counter();ref,mesh,state,parameters,source,time=case_state(case)
        arrays=dict(conserved=state.conserved.copy(),psi=state.psi.copy(),g=state.g.copy(),electron_entropy=state.electron_entropy.copy())
        try:
            rate=rhs(mesh,state,parameters,time=time,source=source)
            metrics,targets=measure(mesh,ref,state,parameters,rate,time,family=case['family'])
        except (RuntimeError,ValueError) as error:
            record=dict(completed=False,failure=repr(error))
        else:
            arrays.update(rate_conserved=rate.state.conserved,rate_psi=rate.state.psi,rate_g=rate.state.g,
                rate_electron_entropy=rate.state.electron_entropy,electric=rate.electric,
                radial_flux=rate.radial_flux,axial_flux=rate.axial_flux,fluid_source=rate.fluid_source,
                density_rate=rate.number_density_rate,regular_velocity=rate.regular_velocity,**targets)
            record=dict(completed=True,failure=None,metrics=metrics,linear_residual=rate.linear_residual,
                linear_iterations=rate.linear_iterations,linear_refinements=rate.linear_refinements)
        for key,array in zip(('conserved','psi','g','electron_entropy'),state.arrays()):np.testing.assert_array_equal(arrays[key],array)
        assert frozen==provenance()
        path=directory/(case['name']+'.npz');np.savez_compressed(path,**arrays)
        record.update(configuration=case,parameters=asdict(parameters),provenance=frozen,elapsed_s=perf_counter()-start,
                      fields=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()),input_unchanged=True)
        write_record(directory/(case['name']+'.json'),record)
        print(case['name'],record.get('metrics'),record['failure'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',default='outputs/finite_mass_initial_attempt');parser.add_argument('--group')
    args=parser.parse_args();run(args.directory,args.group)
