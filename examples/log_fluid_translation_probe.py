"""Quantify the thermal obstruction to an energy-only pressure-force correction."""
import argparse
import json
from pathlib import Path
import numpy as np
from plasma_bench.compatible_bulk_thermal import CompatibleBulkThermal
from .bulk_thermal_verification import initialize
from .log_fluid_records import provenance,sha


def calculate(spec):
    ref,model,flux,_,_,species,states=initialize(spec)
    speed=2000*spec.get('direction',1);momentum=np.zeros((model.mesh.nr,model.mesh.nz,2))
    momentum[...,1]=speed*model.density_map.cell_means
    faces=flux.from_momentum(momentum);nt=flux.density_rate(faces);records={}
    arrays=dict(density=model.density_map.cell_means,momentum=momentum,density_rate=nt,radial=faces.radial,axial=faces.axial)
    for label,thermal,state in zip(('electron','ion'),species,states):
        rate=CompatibleBulkThermal(thermal,flux).rate(state,faces);energy=thermal.energy_rate(state,rate.entropy,nt)
        total=float(np.sum(model.mesh.volumes*energy['total']))
        scale=float(np.sum(model.mesh.volumes*(abs(energy['density'])+abs(energy['entropy']))))
        records[label]=dict(thermal_power_W=total,uncancelled_power_W=scale,relative_chain_rule_defect=abs(total)/scale,
                            force_needed_to_cancel_N=-total/speed)
        arrays.update({label:state,label+'_rate':rate.entropy,label+'_energy_rate':energy['total']})
    return dict(spec=spec,speed_m_per_s=speed,species=records),arrays


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);d=parser.parse_args().directory
    d.mkdir(parents=True,exist_ok=False);records=[]
    identity=provenance();probe=Path(__file__)
    specs=[dict(nr=n,direction=sign) for sign in (1,-1) for n in (8,16,32,64)]+[dict(nr=16,quadrature_order=16)]
    for i,spec in enumerate(specs):
        record,arrays=calculate(spec);path=d/f'translation_{i:02d}.npz';np.savez_compressed(path,**arrays)
        record.update(arrays_file=path.name,arrays_sha256=sha(path));records.append(record);print(spec,record['species'],flush=True)
    (d/'probe_source.py').write_bytes(probe.read_bytes());(d/'probe_source.md').write_bytes(probe.with_suffix('.md').read_bytes())
    (d/'records.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'],probe_sha256=sha(probe),records=records),indent=2)+'\n')


if __name__=='__main__':main()
