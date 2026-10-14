"""Full contact evolution with explicit, resolution-dependent numerical eta6."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME, mu_0 as MU0
from . import guarded_trajectory_controls as guard
from . import matched_mode_evolution as evolution
from .implicit_plasma_records import provenance
from plasma_bench.hyperdiffusive_wall_plasma import HyperdiffusiveWallPlasma
from plasma_bench.hyperdiffusive_modal_gauss import HyperdiffusiveModalGauss


def run(destination, control, spec, *, strength=.25, epsilon=-1e-4, steps=60):
    if not np.isfinite(strength) or strength <= 0:
        raise ValueError('Positive finite numerical strength required')

    def model(system):
        d = system.fluid
        # Inputs are fixed from the original contact, not measured from its error.
        from .continuum_contact_mode import ContinuumContactMode
        ref = ContinuumContactMode().ref
        n0 = ref.number_density
        mass = ref.ion_mass+ME
        sound2 = E*(ref.gamma_e*ref.electron_ev+ref.gamma_i*ref.ion_ev)/mass
        alfven2 = ref.bias**2/(MU0*mass*n0)
        fast = np.sqrt(sound2+alfven2)
        h = min(d.radius/system.geometry.mesh.nr, d.length/system.geometry.mesh.nz)
        delta = (ref.ion_mass-ME)/mass
        hall_diffusivity = abs(delta*ref.bias)/(MU0*E*n0)
        eta6 = MU0*strength*(fast+hall_diffusivity/h)*h**5/np.pi**6
        result = HyperdiffusiveWallPlasma(system.mapping, system.ion_mass,
            eta6=eta6, electron_heat_fraction=1., cache_geometry=system.cache_geometry,
            **system.solver_options)
        result.original_system = system.original_system
        (Path(destination)/'numerical_diffusion.json').write_text(json.dumps(dict(
            strength=strength, eta6_ohm_m5=eta6, h_m=h, fast_speed_m_s=fast,
            hall_diffusivity_m2_s=hall_diffusivity, numerical_electron_heat_fraction=1.,
            formula='mu0*strength*(fast+hall_diffusivity/h)*h^5/pi^6',
            scope='Numerical sixth-order damping, not physical resistivity. h is the smaller mean field spacing; radial cells are nonuniform in r. eta6 vanishes as h^4 in the Hall-dominated limit. Reference is the unmodified collisionless continuum mode.'), indent=2)+'\n')
        return result

    def evolve(destination, spec, **options):
        saved = evolution.RecycledModalGauss, evolution.preconditioned, evolution.provenance
        try:
            evolution.RecycledModalGauss = HyperdiffusiveModalGauss
            evolution.preconditioned = model
            evolution.provenance = lambda seeds: provenance([*seeds, 'examples.hyperdiffusive_contact'])
            return evolution.run(destination, spec, **options)
        finally:
            evolution.RecycledModalGauss, evolution.preconditioned, evolution.provenance = saved

    saved = guard.evolve, guard.provenance
    try:
        guard.evolve = evolve
        guard.provenance = lambda seeds: provenance([*seeds, 'examples.hyperdiffusive_contact'])
        return guard.run(destination, control, spec, epsilon, steps)
    finally:
        guard.evolve, guard.provenance = saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('destination'); parser.add_argument('control')
    parser.add_argument('--strength', type=float, default=.25)
    parser.add_argument('--steps', type=int, default=60)
    parser.add_argument('--elements', type=int, default=8)
    parser.add_argument('--nr', type=int, default=32)
    args = parser.parse_args()
    run(args.destination, args.control,
        dict(elements=args.elements, points=64, field_nr=args.nr, field_nz=64),
        strength=args.strength, steps=args.steps)
