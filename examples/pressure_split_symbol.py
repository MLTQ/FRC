"""Continuum frozen principal-symbol diagnostic for the isolated pressure split."""
from dataclasses import asdict
from pathlib import Path
import hashlib
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from .electron_reference import ElectronReference
from .verification_records import write_record


def symbol(ref, r, z, angle):
    f = ref.fields(r, z); kr, kz = np.cos(angle), np.sin(angle)
    s, k = (r/ref.radius)**2, 2*np.pi/ref.length
    ar = -4*ref.density_contrast*r*(1-s)*np.cos(k*z)/ref.radius**2
    az = -ref.density_contrast*k*(1-s)**2*np.sin(k*z)
    pr, pz = -E*f['density']*r*f['er_over_r'], -E*f['density']*f['ez']
    a = (ar*kz-az*kr)/(E*f['density'])
    c = (kr*(pz-ref.gamma*f['pressure']*az)-kz*(pr-ref.gamma*f['pressure']*ar))/(MU0*E*f['density'])
    b = (f['jr']*kr+f['jz']*kz)/(E*f['density'])
    return a, b, c


def run(directory='outputs/warm_pressure_failure'):
    ref = ElectronReference()
    r = np.linspace(.001, .034, 60); z = np.linspace(-.12, .12, 121); angle = np.linspace(0, 2*np.pi, 361)
    a, b, c = symbol(ref, r[:, None, None], z[None, :, None], angle[None, None, :])
    discriminant = b*b+4*a*c; idx = np.unravel_index(np.argmin(discriminant), discriminant.shape)
    point = [float(r[idx[0]]), float(z[idx[1]]), float(angle[idx[2]])]
    matrix = np.array([[0., a[idx]], [c[idx], b[idx]]])
    eigenvalues = np.linalg.eigvals(1j*matrix)
    growth = .5*np.sqrt(-discriminant[idx])
    np.testing.assert_allclose(max(eigenvalues.real), growth, rtol=1e-14)
    # With uniform n, a=0 and the pressure principal symbol has real speeds.
    uniform = ElectronReference(density_contrast=0.)
    au, bu, cu = symbol(uniform, *point)
    assert au == 0 and bu*bu+4*au*cu >= 0
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    paths = [Path(__file__), Path('examples/electron_reference.py')]
    report = dict(reference=asdict(ref), radial_samples_m=r.tolist(), axial_samples_m=z.tolist(), direction_angles_rad=angle.tolist(),
                  sampled_minimum_discriminant_m2_s2=float(discriminant[idx]), point_r_z_angle=point,
                  real_principal_matrix=matrix.tolist(), growth_rate_per_wavenumber_m_s=float(growth),
                  eigenvalues_i_matrix_real=eigenvalues.real.tolist(), eigenvalues_i_matrix_imag=eigenvalues.imag.tolist(),
                  sources={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  scope='Local frozen pressure-only principal symbol, including background current advection; not a theorem for full warm Hall MHD')
    write_record(directory/'principal_symbol.json', report)
    print(report['sampled_minimum_discriminant_m2_s2'], point, growth, flush=True)


if __name__ == '__main__':
    run()
