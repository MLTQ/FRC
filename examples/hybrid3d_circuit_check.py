"""3D plasma-to-coil reciprocity against the exact mutual inductance of circular filaments."""
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.special import ellipk, ellipe
from scipy.constants import mu_0 as MU0
from plasma_bench.hybrid3d import curl_edges
from plasma_bench.hybrid3d_device import Device, device_box, cell_types, masks
from examples.hybrid3d_circuit import CircuitDrive
from plasma_bench.hybrid3d_coils import VacuumCoilBasis


def loop_edges(box, radius, current):
    """Biot–Savart vector potential, independent of the production coil mapper."""
    result = []
    for a in range(2):
        axes = [box.lo[q]+box.d[q]*(np.arange(box.shape[q])+(.5 if q == a else 0.)) for q in range(3)]
        x, y, z = np.meshgrid(*axes, indexing='ij')
        r = np.hypot(x, y)
        m = np.clip(4*radius*r/((radius+r)**2+z*z), 1e-20, 1-1e-13)
        ap = MU0*current/np.pi*np.sqrt(radius/np.maximum(r*m, 1e-30))*((1-m/2)*ellipk(m)-ellipe(m))
        ap = np.where(r < 1e-10, 0., ap)
        result.append(ap*(-y if a == 0 else x)/np.maximum(r, 1e-30))
    return result+[np.zeros(box.shape)]


def reciprocity(dx, inner=0., project_source=False):
    dev = Device(.14, .2, .18, .03, inner)
    box = device_box(dev, dx, dx)
    ap, ac, current = .08137, .16537, 1e4
    source = np.array(loop_edges(box, ac, 1.))
    if project_source:
        source = VacuumCoilBasis(box, cell_types(box, dev) == 0)(source)
    plasma = loop_edges(box, ap, current)
    b = list(curl_edges(*plasma, box.d))
    model = SimpleNamespace(b=b, d=box.d, box=box, edge_mask=masks(cell_types(box, dev))[2], time=0.)
    drive = CircuitDrive.__new__(CircuitDrive)
    drive.edges, drive.faces = source[None], np.zeros((1, 3)+box.shape)
    drive.coefficients = lambda t: np.zeros(1)
    measured = float(drive.linkage(model)[0])
    m = 4*ap*ac/(ap+ac)**2
    mutual = MU0*np.sqrt(ap*ac)/np.sqrt(m)*((2-m)*ellipk(m)-2*ellipe(m))
    exact = current*mutual
    return dict(dx_m=dx, inner_m=inner, projected_source=project_source, exact_linkage_wb=exact, measured_linkage_wb=measured,
                relative_error=(measured-exact)/exact)


def main(out='outputs/hybrid3d_circuit_check'):
    out = Path(out);out.mkdir(exist_ok=True, parents=True)
    result = [reciprocity(dx, inner) for inner in (0., .04) for dx in (.02, .01, .005)]
    (out/'reciprocity.json').write_text(json.dumps(result, indent=1))
    print(json.dumps(result, indent=1), flush=True)
    return result


if __name__ == '__main__':
    main()
