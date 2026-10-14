"""Explicit tube, pre-ionized inventory and biased-waveform material experiment."""
import numpy as np
from scipy.constants import atomic_mass, Boltzmann as KB, torr
from .mhd_mesh import MHDMesh
from .material_cylinder import MaterialCylinder
from .material_two_species import MaterialTwoSpecies
from .material_vacuum_space import MaterialVacuumSpace
from .material_magnetic_vacuum import MaterialMagneticVacuum
from .material_coil_current import MaterialCoilCurrent
from .electromagnetics import loop_field
from .waveforms import from_dict


# Ion masses are declared natural-abundance approximations. H2 contributes two
# protons per initial neutral molecule after the assumed complete dissociation.
GASES = {'Ar': (39.948, 1), 'He': (4.002602, 1),
         'Ne': (20.1797, 1), 'H2': (1.00794, 2)}


class BiasedWaveform:
    """Common coil current corresponding to a DC field plus a serialized pulse."""
    def __init__(self, bias, pulse, field_per_amp):
        self.bias = bias
        self.pulse = pulse
        self.field_per_amp = field_per_amp

    def value(self, time):
        return (self.bias + self.pulse.value(time)) / self.field_per_amp

    def derivative(self, time):
        return self.pulse.derivative(time) / self.field_per_amp


class TubeMaterialCoilCurrent(MaterialCoilCurrent):
    """Free evolution only until contact with the explicitly nonconducting tube."""
    def __init__(self, *args, tube_radius, tube_half_length, **kwargs):
        super().__init__(*args, **kwargs)
        self.tube_radius = tube_radius
        self.tube_half_length = tube_half_length

    def response(self, state, *, rate=True, time=0.):
        boundary = state.fluid.coordinates[self.fluid.domain.boundary]
        if np.any(boundary[:, 0] >= self.tube_radius) or np.any(abs(boundary[:, 1]) >= self.tube_half_length):
            raise RuntimeError('Plasma reached tube: wall/sheath/contact physics is not implemented')
        return super().response(state, rate=rate, time=time)


def build(config):
    """Build the actual solver; fill pressure fixes inventory, not local density."""
    tube = config['tube']; plasma = config['plasma']; grid = config['mesh']
    source = config['source']; evolution = config['evolution']
    values = [tube['radius_m'], tube['length_m'], plasma['radius_m'],
              plasma['half_length_m'], plasma['neutral_fill_mtorr'],
              plasma['neutral_temperature_k'], *plasma['central_temperature_ev'],
              evolution['duration_s']]
    if not np.isfinite(values).all() or min(values) <= 0 or len(plasma['central_temperature_ev']) != 2:
        raise ValueError('Positive finite dimensions, fill, temperatures and duration required')
    if plasma['assumed_charge_state'] != 1 or plasma['assumed_ionization_fraction'] != 1.:
        raise ValueError('Only an explicitly assumed fully singly ionized inventory is implemented')
    if plasma['radius_m'] >= tube['radius_m'] or plasma['half_length_m'] >= tube['length_m']/2:
        raise ValueError('Initial free plasma must fit strictly inside the nonconducting tube')
    if not isinstance(evolution['steps'], int) or evolution['steps'] < 4:
        raise ValueError('At least four time steps required')
    if not np.isfinite(source['bias_center_t']):
        raise ValueError('Finite prescribed bias required')
    mass, ions_per_neutral = GASES[plasma['gas']]
    mesh = MHDMesh(plasma['radius_m'], plasma['half_length_m'], grid['nr'], grid['nz'])
    domain = MaterialCylinder(mesh, quadrature=grid['quadrature'])
    space = MaterialVacuumSpace(mesh, degree=grid['degree'], quadrature=grid['quadrature'])
    r, z = space.reference_r, space.reference_z
    profile = ((1-(r/mesh.radius)**2)*(1-(z/mesh.half_length)**2)).reshape(mesh.nr*mesh.nz, -1)
    cell_area = (np.diff(mesh.rf)[:, None]*np.diff(mesh.zf)[None, :]).ravel()
    volume = 2*np.pi*r.reshape(profile.shape)*cell_area[:, None]*domain.qweight[None, :]
    chamber_volume = np.pi*tube['radius_m']**2*tube['length_m']
    neutral_density = plasma['neutral_fill_mtorr']*1e-3*torr/(KB*plasma['neutral_temperature_k'])
    inventory = ions_per_neutral*neutral_density*chamber_volume
    center_density = inventory/np.sum(volume*profile**1.5)
    temperatures = np.asarray(plasma['central_temperature_ev'])
    fluid = MaterialTwoSpecies(domain, domain.fluid_reference, center_density*profile**1.5,
                               profile[..., None]*temperatures, mass*atomic_mass)
    # The profile is isentropic initially. Set the analytically constant values
    # directly to avoid floating roundoff creating unsupported dual-volume modes.
    fluid.initial.entropy[:] = np.log(temperatures)-(2/3)*np.log(center_density)
    vacuum = MaterialMagneticVacuum(mesh, degree=grid['degree'], quadrature_order=grid['quadrature'],
                                    vacuum_cells=grid['vacuum_cells'], modes=grid['vacuum_modes'])
    coils = tuple(tuple(row) for row in source['coils_radius_z_turns'])
    for radius, center, turns in coils:
        if not np.isfinite([radius, center, turns]).all() or radius <= tube['radius_m']:
            raise ValueError('Finite prescribed loops must lie radially outside the tube')
    field_per_amp = sum(turns*float(loop_field(0., 0., radius, center)[2]) for radius, center, turns in coils)
    if not np.isfinite(field_per_amp) or field_per_amp == 0:
        raise ValueError('Nonzero finite vacuum-center field per common ampere required')
    pulse = from_dict(source['pulse_center_t'])
    if float(pulse.value(0.)) != 0.:
        raise ValueError('Pulse must start at zero so the initial field equals the declared DC bias')
    wave = BiasedWaveform(source['bias_center_t'], pulse, field_per_amp)
    model = TubeMaterialCoilCurrent(fluid, space, vacuum, coils, wave,
        resistivity=plasma['prescribed_resistivity_ohm_m'], coulomb_exchange=True,
        tube_radius=tube['radius_m'], tube_half_length=tube['length_m']/2)
    metadata = dict(central_density_m3=float(center_density), ion_inventory=float(fluid.particles.sum()),
        assumed_neutral_inventory=float(neutral_density*chamber_volume), ions_per_neutral=ions_per_neutral,
        neutral_density_before_preparation_m3=float(neutral_density), ion_mass_kg=float(fluid.ion_mass),
        vacuum_center_field_per_amp_t=float(field_per_amp), actual_initial_entropy=fluid.initial.entropy[0, 0].tolist(),
        assumptions=['All initial gas inventory is assumed prepared as the declared compact singly ionized plasma.',
            'No ionization, recombination, neutrals, radiation or preparation energy is evolved.',
            'Quartz is represented only as a nonconducting, nonmagnetic contact limit; no wall or sheath model.',
            'Scalar resistivity is prescribed. Coulomb temperature exchange remains conditional.',
            'Cold dilute edge, first-order entropy diffusion and current time sensitivity remain unresolved.'])
    return fluid, model, metadata
