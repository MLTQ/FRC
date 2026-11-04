"""Local experiment definitions, validated launch inputs and existing FRC adapter."""
import json
from pathlib import Path
import math

from examples.iec_3d import DEFAULTS


PARAMETERS = [
    ('bias_v', 'Cage bias', 'V', 0., 3000., 10.),
    ('duration_ns', 'Duration', 'ns', 1., 2000., 10.),
    ('particles', 'Macroparticles / species', '', 100, 10000, 100),
    ('source_ua', 'Wall electron source', 'µA', 0., 1000., 10.),
    ('source_ev', 'Emitted electron energy', 'eV', .01, 100., .1),
    ('density_m3', 'Initial seed density', 'm⁻³', 1e8, 1e14, 1e10),
    ('electron_ev', 'Initial electron temperature', 'eV', .01, 100., .1),
    ('ion_ev', 'Initial ion temperature', 'eV', .001, 10., .01),
    ('ion_amu', 'Ion mass (hydrogenic, Z=1)', 'mₚ', 1., 4., 1.),
    ('rise_ns', 'Bias rise time', 'ns', 0., 1000., 1.),
    ('modulation', 'Sinusoidal modulation depth', '', 0., 1., .05),
    ('frequency_mhz', 'Modulation frequency', 'MHz', .001, 100., .1),
    ('dt_ps', 'Integration timestep', 'ps', 1., 100., 1.),
    ('mesh_mm', 'Bulk mesh target', 'mm', 8., 30., 1.),
    ('wire_mesh_mm', 'Wire mesh target', 'mm', .5, 1.5, .1),
]


def catalog():
    return [dict(id='iec', name='Positive cage', family='Electrostatic · kinetic',
                 description='Full 3D electrons + ions, self-consistent space charge and absorbing steel.',
                 scope='Collisionless development model. Spherical chamber; idealized wall electron emitter. No air chemistry, secondary emission, reactions or fusion gain.',
                 parameters=[dict(key=k, label=l, unit=u, min=a, max=b, step=s, default=DEFAULTS[k]) for k,l,u,a,b,s in PARAMETERS]),
            dict(id='frc', name='FRC formation', family='Magnetic · hybrid / MHD',
                 description='Existing formation runs, full 3D hybrid particles, magnetic fields and current density.',
                 scope='Existing scientific viewer and source-frozen campaigns. Physical edge losses and 3D convergence remain unresolved.', parameters=[])]


def validate(values):
    if not isinstance(values, dict) or set(values)-{p[0] for p in PARAMETERS}:
        raise ValueError('Unknown experiment parameter')
    config = DEFAULTS.copy()
    for key, label, _, minimum, maximum, _ in PARAMETERS:
        value = values.get(key, config[key])
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not minimum <= value <= maximum:
            raise ValueError(f'{label} must be between {minimum:g} and {maximum:g}')
        if key == 'particles':
            if int(value) != value:
                raise ValueError('Particle count must be an integer')
            value = int(value)
        config[key] = value
    if config['duration_ns']*1000/config['dt_ps'] > 200000:
        raise ValueError('This local launcher supports at most 200,000 timesteps per run')
    return config


def runs(root):
    root = Path(root)
    rows = []
    paths = list((root/'outputs/lab_runs').glob('*/run.json'))+list((root/'outputs').glob('iec_3d_*/run.json'))
    for path in sorted(paths, reverse=True):
        try:
            row = json.loads(path.read_text())
            row['url'] = '/data/'+str(path.parent.relative_to(root/'outputs'))
            rows.append(row)
        except (OSError, ValueError):
            continue
    source = root/'outputs/formation_display/runs.json'
    if source.exists():
        for row in json.loads(source.read_text()):
            rows.append(dict(id=row['run'], experiment='frc', label=row['label'],
                             status='saved', frames=row.get('frames', 0),
                             description=row.get('description', ''),
                             viewer_url='/data/formation_display/viewer.html?run='+row['run']))
    return rows
