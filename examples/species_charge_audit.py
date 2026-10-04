"""Replay the pressure/entropy probe and check its continuum work identity."""
import hashlib
import json
from pathlib import Path
from .species_charge_probe import measure
import numpy as np


def main():
    root=Path(__file__).resolve().parents[1];out=root/'outputs/species_charge_probe';source=root/'outputs/current_transport_verification'
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    provenance=json.loads((out/'provenance.json').read_text());summary=json.loads((out/'summary.json').read_text())
    for name,digest in provenance['source_files'].items():assert sha(root/name)==digest,name
    results=[]
    for record in summary['records']:
        raw=source/Path(record['base_record']).with_suffix('.npz');assert sha(raw)==record['arrays_sha256'];arrays=dict(np.load(raw))
        for order in (16,20):assert measure(record['spec'],arrays,order)==record[f'order{order}']
        checks=[]
        for lo,hi in zip(record['order16'],record['order20']):
            scale=hi['work_scale'];identity=abs(hi['unpaired_pressure_work']-hi['physical_gradient_work'])/scale
            compatible=abs(hi['compatible_gradient_work'])/scale;quadrature=abs(hi['pressure_work']-lo['pressure_work'])/scale
            assert max(identity,compatible,quadrature)<1e-13
            checks.append(dict(species=hi['species'],work_identity_scaled=identity,compatible_gradient_scaled=compatible,pressure_quadrature_scaled=quadrature))
        results.append(dict(spec=record['spec'],orders16_and20_bitwise=True,checks=checks));print(record['spec'],'audited',flush=True)
    for k in range(2):
        values=[r['order20'][k]['potential_reconstruction_relative'] for r in summary['records']]
        assert all(a>b for a,b in zip(values,values[1:]))
    result=dict(source_sha256=provenance['source_sha256'],auditor_sha256=sha(Path(__file__)),records=results)
    (out/'audit.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
