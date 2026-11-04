"""Actual open-cage spatial refinement; cross-check independent vacuum solver."""
import argparse
import json
from pathlib import Path

import numpy as np

from plasma_bench.iec_field import Field
from plasma_bench.iec_mesh import build_mesh


def run(out, reference=None):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for size, wire in ((.024, .0015), (.018, .001), (.014, .00075)):
        path = out/f'cage_{size}_{wire}.npz'
        if reference and size == .018:
            import shutil
            shutil.copy2(reference, path)
            mesh = dict(np.load(path))
        else:
            mesh = build_mesh(path, size=size, wire_size=wire)
        field = Field(mesh)
        centre = np.zeros((1, 3))
        cell = field.locate(centre)
        phi0 = float(field.weights(centre, cell)[0]@field.unit[field.tets[cell[0]]])
        row = dict(bulk_size_m=size, wire_size_m=wire, nodes=len(field.nodes), tetrahedra=len(field.tets),
                   capacitance_pf=field.capacitance*1e12, centre_phi_per_bias=phi0,
                   boundary_reference_capacitance_pf=11.090, boundary_reference_phi0=.8660754,
                   note='Independent MFS reference has up to 1.4% local surface potential error; neither solution is exact.')
        rows.append(row)
        (out/'report.json').write_text(json.dumps(rows, indent=2, allow_nan=False))
        print(row, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output')
    parser.add_argument('--reference', type=Path)
    args = parser.parse_args()
    run(args.output, args.reference)
