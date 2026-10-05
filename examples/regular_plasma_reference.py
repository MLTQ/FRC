"""Independent continuum wave/contact initialization on the actual regular state."""
import numpy as np
from scipy.constants import atomic_mass, electron_mass, elementary_charge as E
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.regular_field_map import RegularFieldMap
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from plasma_bench.regular_plasma import RegularPlasma
from .finite_mass_contact_reference import FiniteMassContact
from .finite_mass_wave_reference import FiniteMassWave
from .finite_mass_manufactured import cell_average


def reference(spec):
    if spec['kind'] == 'contact':
        return FiniteMassContact(bias=spec.get('bias', .005))
    if spec['kind'] != 'wave':
        raise ValueError('Known independent contact or wave required')
    n = spec.get('density', 1e18)
    return FiniteMassWave(radius=.035, length=.24, number_density=n,
        ion_mass=39.948*atomic_mass, bias=spec.get('bias', .005), pressure=6*n*E,
        gamma=5/3, gamma_i=1.4, branch=spec.get('branch', 2),
        amplitude=spec.get('amplitude', .1))


def exact(ref, kind, r, z, time=0.):
    """Physical variables and independently differentiated time rates."""
    r, z = np.broadcast_arrays(r, z)
    if kind == 'contact':
        f = ref.fields(z, time)
        t = ref.fields(z, time, rate=True)
        n, nt = f['number_density'], t['number_density']
        v = np.zeros(r.shape+(3,)); v[..., 2] = ref.velocity
        electric = np.zeros_like(v); electric[..., 2] = f['electric_z']
        vt = np.zeros_like(v); b = np.zeros_like(v); bt = b.copy()
        current = np.zeros_like(v); jt = current.copy()
    else:
        f = ref.fields(r, z, time)
        t = {k: v.real for k, v in ref.derivatives(r, z, time)[2].items()}
        n = ref.number_density+f['density']/(ref.ion_mass+ref.electron_mass)
        nt = t['density']/(ref.ion_mass+ref.electron_mass)
        v, vt = f['velocity'], t['velocity']
        electric, b, bt = f['electric'], f['magnetic'], t['magnetic']
        current, jt = f['current'], t['current']
    p = np.stack((f['electron_pressure'], f['ion_pressure']), axis=-1)
    pt = np.stack((t['electron_pressure'], t['ion_pressure']), axis=-1)
    if kind == 'wave':
        p += np.array([ref.electron_pressure, ref.ion_pressure])
    return dict(n=n, nt=nt, velocity=v, velocity_rate=vt, pressure=p, pressure_rate=pt,
                electric=electric, magnetic=b, magnetic_rate=bt, current=current, current_rate=jt)


def fluid_state(ref, kind, fluid, time=0.):
    f = exact(ref, kind, fluid.r, fluid.z, time)
    n = f['n']; v = f['velocity']
    u = np.zeros(fluid.shape)
    u[..., 0] = n
    np.divide(n[..., None]*v[..., :2], fluid.r[..., None], out=u[..., 1:3],
              where=fluid.r[..., None] > 0)
    u[..., 3] = n*v[..., 2]
    s = np.log(f['pressure']/(E*n[..., None]))
    s -= (fluid.thermal.gammas-1)*np.log(n/fluid.thermal.n0)[..., None]
    u[..., 4:] = n[..., None]*s
    return u


def initialize(spec):
    ref = reference(spec); ne = spec['elements']
    gammas = (ref.gamma_e, ref.gamma_i) if spec['kind'] == 'contact' else (ref.gamma, ref.gamma_i)
    d = RegularCoordinateFluid(ref.radius, ref.length, ne, spec.get('points', 8*ne),
        ref.ion_mass+electron_mass, degree=spec.get('degree', 5), gammas=gammas)
    mesh = MHDMesh(ref.radius, .5*ref.length, spec.get('field_nr', max(8, 2*ne)),
                   spec.get('field_nz', max(8, 2*ne)))
    geometry = CompatibleGeometry(mesh, degree=5, quadrature_order=spec.get('quadrature', 12))
    mapping = RegularFieldMap(d, geometry)
    system = RegularPlasma(mapping, ref.ion_mass, tolerance=spec.get('tolerance', 1e-10))
    u = fluid_state(ref, spec['kind'], d)
    psi = np.zeros((mesh.nr+1, mesh.nz+1)); g = np.zeros((mesh.nr, mesh.nz))
    if spec['kind'] == 'wave':
        psi = ref.fields(*mesh.corners)['psi']; psi[[0, -1]] = 0.; psi[:, -1] = psi[:, 0]
        g = cell_average(mesh, lambda r, z: ref.fields(r, z)['magnetic'][..., 1]/r, order=16)
    field = HarmonicMagnetic(mesh, psi, g, ref.bias)
    return ref, d, u, system, field
