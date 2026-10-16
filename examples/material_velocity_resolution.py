"""Assess velocity coarsening on actual failed reacting-material states."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.constants import atomic_mass, electron_mass as ME
from scipy.linalg import eigh
from .implicit_plasma_records import provenance, freeze, sha


def restrict_velocity(faces, numbers, mass, velocity, thermal, count):
    """Nested continuous P1 velocity projection; density/thermal cuts stay put."""
    target=faces[0]+(faces[-1]-faces[0])*(1-np.cos(np.linspace(0,np.pi,count+1)))/2
    selected=np.unique(np.argmin(abs(faces[:,None]-target[None,:]),axis=0))
    coarse=faces[selected];h=np.diff(faces);cell=np.searchsorted(coarse,faces[:-1],side='right')-1
    width=np.diff(coarse)[cell];a=(faces[:-1]-coarse[cell])/width;b=(faces[1:]-coarse[cell])/width
    matrix=np.zeros((len(coarse),len(coarse)));stiffness=np.zeros_like(matrix);load=np.zeros((len(coarse),3))
    cell_mass=mass*numbers;second=(a*a+a*b+b*b)/3;first=(a+b)/2
    for j,m,s1,s2,k in zip(cell,cell_mass,first,second,(5/3)*(2/3)*thermal/width**2):
        matrix[j:j+2,j:j+2]+=m*np.array([[1-2*s1+s2,s1-s2],[s1-s2,s2]])
        stiffness[j:j+2,j:j+2]+=k*np.array([[1.,-1.],[-1.,1.]])
    vl,vr=velocity[:-1],velocity[1:]
    np.add.at(load,cell,cell_mass[:,None]*((3-2*a-b)[:,None]*vl+(3-a-2*b)[:,None]*vr)/6)
    np.add.at(load,cell+1,cell_mass[:,None]*((2*a+b)[:,None]*vl+(a+2*b)[:,None]*vr)/6)
    scale=np.max(np.diag(matrix));projected=np.linalg.solve(matrix/scale,load/scale)
    left=(1-a[:,None])*projected[cell]+a[:,None]*projected[cell+1]
    right=(1-b[:,None])*projected[cell]+b[:,None]*projected[cell+1]
    kinetic=lambda l,r:cell_mass/6*np.sum(l*l+l*r+r*r,axis=1)
    old_kinetic=kinetic(vl,vr);new_kinetic=kinetic(left,right);recovered=thermal+old_kinetic-new_kinetic
    momentum_change=np.sum(cell_mass[:,None]*(left+right-vl-vr)/2,axis=0)
    try:frequency=float(np.sqrt(eigh(stiffness/scale,matrix/scale,subset_by_index=[len(coarse)-1]*2,eigvals_only=True)[0]))
    except np.linalg.LinAlgError:frequency=None
    negative=recovered<0
    record=dict(density_cells=len(numbers),velocity_cells=len(coarse)-1,minimum_density_width_m=float(h.min()),
        minimum_velocity_width_m=float(np.diff(coarse).min()),minimum_recovered_thermal_j_m2=float(recovered.min()),
        inadmissible_cells=int(negative.sum()),inadmissible_number_fraction=float(numbers[negative].sum()/numbers.sum()),
        negative_thermal_sum_j_m2=float(recovered[negative].sum()),
        total_unresolved_kinetic_j_m2=float((old_kinetic-new_kinetic).sum()),
        momentum_change_kg_per_m_s=momentum_change.tolist(),
        acoustic_angular_frequency_bound_s=frequency,
        maximum_recovered_thermal_deficit_fraction=float(np.max(np.divide(-recovered,thermal+old_kinetic,out=np.zeros_like(thermal),where=(thermal+old_kinetic)>0),initial=0.)),
        admissible=not bool(negative.any()))
    arrays=dict(density_faces_m=faces,numbers_m2=numbers,original_velocity_m_s=velocity,original_thermal_j_m2=thermal,
        coarse_faces_m=coarse,projected_velocity_m_s=projected,recovered_thermal_j_m2=recovered,mass_matrix=matrix,acoustic_stiffness=stiffness)
    return record,arrays


def run(source_directory,destination):
    source_directory=Path(source_directory);out=Path(destination)
    if out.exists():raise ValueError('Fresh velocity-resolution output required')
    source=json.loads((source_directory/'summary.json').read_text())
    identity=provenance(['examples.material_velocity_resolution']);out.mkdir(parents=True);freeze(out,identity)
    records=[];mass=np.r_[39.948*atomic_mass-np.arange(19)*ME,ME]
    for case in source['cases']:
        count=int(case['name'].split('n')[-1]);raw=source_directory/(case['name']+'.npz');values=np.load(raw);step=case['steps']
        for block,species in enumerate(values[f'step{step}_species_indices']):
            prefix=f'step{step}_block{block}_';record,arrays=restrict_velocity(values[prefix+'faces_m'][0],values[prefix+'numbers_m2'][0],mass[species],
                values[prefix+'relative_velocity_m_s'][0],values[prefix+'thermal_j_m2'][0],count)
            record.update(case=case['name'],species=int(species),block=block,time_s=case['completed_time_s'],input_sha256=sha(raw))
            np.savez_compressed(out/f'{case["name"]}_species{species}_block{block}.npz',**arrays);records.append(record)
    report=dict(source_sha256=identity['source_sha256'],input_source_sha256=source['source_sha256'],records=records,
        scope='Snapshot assessment of a proposed coarser velocity space. Actual density cuts, species support and charge are unchanged. No state is advanced or adopted; negative recovered thermal energy rejects the candidate. Acoustic spectrum uses the original adiabatic pressure, excludes electric/collision sources, and is not a full-system timestep certificate.')
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    for case in source['cases']:
        rows=[r for r in records if r['case']==case['name']]
        print(case['name'],dict(blocks=len(rows),inadmissible_blocks=sum(not r['admissible'] for r in rows),
            inadmissible_cells=sum(r['inadmissible_cells'] for r in rows),maximum_deficit_fraction=max(r['maximum_recovered_thermal_deficit_fraction'] for r in rows)),flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source_directory');parser.add_argument('destination');args=parser.parse_args()
    run(args.source_directory,args.destination)
