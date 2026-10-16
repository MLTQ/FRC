"""Spatial species transport with longitudinal Gauss and transverse Maxwell fields."""
import numpy as np
from scipy.constants import c as LIGHT_SPEED, elementary_charge as EV, epsilon_0 as EPS0
from .atomic_electrostatic import AtomicElectrostatic


class AtomicMaxwell(AtomicElectrostatic):
    """Reflecting material walls and perfectly conducting electromagnetic ends.

    Transverse E lives at faces and transverse B at cell centres. Bz is a
    spatial constant; the longitudinal electric field follows the actual charge.
    All spatial coordinates and fields are in the laboratory frame.
    """
    def __init__(self, *args, axial_magnetic=0., **kwargs):
        super().__init__(*args, **kwargs)
        if not np.isfinite(axial_magnetic):
            raise ValueError('Finite constant axial magnetic field required')
        self.axial_magnetic = float(axial_magnetic)

    def validate_fields(self, electric, magnetic):
        if (electric.shape != (self.cells+1, 2) or magnetic.shape != (self.cells, 2)
                or not np.isfinite(electric).all() or not np.isfinite(magnetic).all()):
            raise ValueError('Finite face transverse E and cell transverse B required')
        if np.any(electric[[0, -1]] != 0):
            raise ValueError('Conducting ends require zero tangential E, without projection')

    def magnetic_vector(self, transverse):
        return np.column_stack((transverse, np.full(self.cells, self.axial_magnetic)))

    def number_velocity(self, state):
        """Return f_s v_s in laboratory units, also for exactly empty species."""
        s = self.source
        return (s.momentum_scale/s.mass[:, None])*state[:, s.momentum] + (
            state[:, self.number, None]*self.frame*s.speed_scale)

    def timestep(self, state, electric, magnetic):
        self.validate_fields(electric, magnetic)
        dt = min(super().timestep(state), self.cfl*self.dx/LIGHT_SPEED)
        gyro = EV*abs(self.charge)/self.source.mass
        occupied = state[:, self.number] > 0
        omega = np.max(np.where(occupied, gyro*np.linalg.norm(self.magnetic_vector(magnetic), axis=1)[:, None], 0.))
        return min(dt, self.maximum_plasma_phase/omega) if omega else dt

    def coupled_derivative(self, state, electric, magnetic):
        self.validate_fields(electric, magnetic)
        rate = super().derivative(state)
        s = self.source
        centre_e = .5*(electric[:-1]+electric[1:])
        e_vector = np.column_stack((centre_e, np.zeros(self.cells)))
        nv = self.number_velocity(state)
        force = (EV/s.momentum_scale)*self.charge[:, None]*(
            state[:, self.number, None]*e_vector[:, None, :]
            + np.cross(nv, self.magnetic_vector(magnetic)[:, None, :]))
        rate[:, s.momentum] += force
        work = self.charge*np.sum(nv[:, :, :2]*centre_e[:, None, :], axis=2)
        rate[:, s.energy] += work - (s.momentum_scale/EV)*s.speed_scale*np.sum(force*self.frame, axis=2)

        cell_current = EV*self.density*np.einsum('s,csi->ci', self.charge, nv[:, :, :2])
        face_current = .5*(cell_current[:-1]+cell_current[1:])
        gradient_b = np.diff(magnetic, axis=0)/self.dx
        gradient_e = np.diff(electric, axis=0)/self.dx
        electric_rate = np.zeros_like(electric)
        electric_rate[1:-1] = LIGHT_SPEED**2*np.column_stack((-gradient_b[:, 1], gradient_b[:, 0])) - face_current/EPS0
        magnetic_rate = np.column_stack((gradient_e[:, 1], -gradient_e[:, 0]))

        field_power = EPS0/(EV*self.density)*self.dx*(
            np.sum(electric[1:-1]*electric_rate[1:-1])
            + LIGHT_SPEED**2*np.sum(magnetic*magnetic_rate))
        self.last_diagnostic['transverse_work_residual_ev_m_s'] = float(self.dx*work.sum()+field_power)
        return rate, electric_rate, magnetic_rate

    def step(self, state, electric, magnetic, time, dt, *, reactions=True, heat=None):
        """Advance the complete coupled spatial system by RK4 and source halves."""
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('Positive finite coupled step required')

        def material(values, start, duration):
            initial = np.zeros(self.cells) if heat is None else np.broadcast_to(heat(self.x, start), (self.cells,))
            middle = np.zeros(self.cells) if heat is None else np.broadcast_to(heat(self.x, start+duration/2), (self.cells,))
            return np.asarray([self.source.step(row, duration, self.density, initial_heat=initial[j], midpoint_heat=middle[j])
                               for j, row in enumerate(values)])

        first = material(state, time, dt/2) if reactions else state
        values = (first, electric, magnetic)
        k1 = self.coupled_derivative(*values)
        k2 = self.coupled_derivative(*(u+.5*dt*k for u, k in zip(values, k1)))
        k3 = self.coupled_derivative(*(u+.5*dt*k for u, k in zip(values, k2)))
        k4 = self.coupled_derivative(*(u+dt*k for u, k in zip(values, k3)))
        result = tuple(u+dt/6*(a+2*b+2*d+e) for u, a, b, d, e in zip(values, k1, k2, k3, k4))
        self.primitive(result[0]); self.validate_fields(*result[1:])
        if reactions:
            result = (material(result[0], time+dt/2, dt/2), *result[1:])
        return result

    def totals(self, state, electric, magnetic):
        self.validate_fields(electric, magnetic)
        row = super().totals(state)
        scale = EPS0/(2*EV*self.density)
        electric_energy = scale*np.trapezoid(np.sum(electric**2, axis=1), dx=self.dx)
        magnetic_energy = scale*LIGHT_SPEED**2*self.dx*np.sum(magnetic**2)
        exact_linear_e = scale*self.dx/3*np.sum(electric[:-1]**2+electric[:-1]*electric[1:]+electric[1:]**2)
        row.update(transverse_electric_energy_ev=electric_energy, transverse_magnetic_energy_ev=magnetic_energy,
                   total_energy_ev=row['total_energy_ev']+electric_energy+magnetic_energy,
                   linear_field_total_energy_ev=row['linear_field_total_energy_ev']+exact_linear_e+magnetic_energy)
        # Constant Bz energy is time independent and excluded from the wave budget.
        return row
