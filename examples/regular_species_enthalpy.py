"""Independent continuous local thermal/current-enthalpy balance audit."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as E
from .regular_species_reference import initialize
from .current_transport_sampling import TransportSampler
from .radial_work_records import freeze, sha


def measure(spec, arrays, order):
    ref, d, _, m, _ = initialize(spec)
    s = TransportSampler(d, m.mesh, order)
    u = arrays['state']
    n = s.scalar(u[..., 0])
    j, dj = s.electric(arrays['current'], True)
    divj = np.trace(dj, axis1=-2, axis2=-1)
    result = []
    for k, (gamma, a) in enumerate(zip(d.thermal.gammas, (-m.beta, m.alpha))):
        entropy = s.scalar(u[..., k+4])/n
        h = E*np.exp(entropy)*(n/d.thermal.n0)**(gamma-1)/(gamma-1)
        dh = []
        for axis in (0, 2):
            dn = s.scalar(u[..., 0], axis=axis)
            ds = (s.scalar(u[..., k+4], axis=axis)-entropy*dn)/n
            dh.append(h*((gamma-1)*dn/n+ds))
        thermal = h*s.scalar(arrays['fluid_rate'][..., k+4])
        electrical = np.sum(j*s.electric(arrays['pressure_electric'][k]), axis=-1)
        # Original current enthalpy flux is a*gamma*h*J/e, including its actual
        # strong divergence. No numerical continuity defect is silently removed.
        divergence = a*gamma/E*(j[..., 0]*dh[0]+j[..., 2]*dh[1]+h*divj)
        defect = thermal-electrical+divergence
        scale2 = thermal**2+electrical**2+divergence**2
        norm = lambda mask: float(np.sqrt(np.sum(s.weight*mask*defect**2)/np.sum(s.weight*mask*scale2)))
        rms = np.sqrt(np.sum(s.weight*scale2)/np.sum(s.weight))
        result.append(dict(local_balance_relative_L2=norm(1.),
                           axis_element=norm(s.R < d.radius/d.elements),
                           wall_element=norm(s.R > d.radius*(1-1/d.elements)),
                           maximum_over_term_rms=float(np.max(abs(defect))/rms),
                           integrated_thermal_rate=float(np.sum(s.weight*thermal)),
                           integrated_electrical_work=float(np.sum(s.weight*electrical)),
                           integrated_enthalpy_divergence=float(np.sum(s.weight*divergence)),
                           integrated_balance_defect=float(np.sum(s.weight*defect))))
    return result


def main(source, output):
    source, output = Path(source), Path(output)
    root = Path(__file__).resolve().parents[1]
    identity = json.loads((source/'provenance.json').read_text())
    files = dict(identity['source_files'])
    files['examples/regular_species_enthalpy.py'] = sha(__file__)
    for name, expected in files.items():
        assert sha(root/name) == expected, name
    audit_identity = dict(identity, source_files=files,
                          source_sha256=hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest(),
                          calculation_source_sha256=identity['source_sha256'])
    output.mkdir(parents=True, exist_ok=False)
    freeze(output, audit_identity)
    records = []
    for i in range(14):
        path = source/f'case_{i:02d}.json'
        row = json.loads(path.read_text())
        raw = source/row['arrays_file']
        assert row['source_sha256'] == identity['source_sha256'] and sha(raw) == row['arrays_sha256']
        arrays = dict(np.load(raw))
        result = dict(file=path.name, spec=row['spec'], arrays_sha256=row['arrays_sha256'],
                      order16=measure(row['spec'], arrays, 16), order20=measure(row['spec'], arrays, 20))
        records.append(result)
        (output/path.name).write_text(json.dumps(result, indent=2)+'\n')
        print(path.name, [r['local_balance_relative_L2'] for r in result['order20']], flush=True)
    for name, expected in files.items():
        assert sha(root/name) == expected, name
    (output/'summary.json').write_text(json.dumps(dict(source_sha256=audit_identity['source_sha256'],
        calculation_source_sha256=identity['source_sha256'], terminal=True, records=records), indent=2)+'\n')


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('source')
    p.add_argument('output')
    args = p.parse_args()
    main(args.source, args.output)
