"""Uniform-background block inverse used only to accelerate a nonlinear solve."""
import numpy as np
from scipy.constants import elementary_charge as E
from scipy.sparse import diags, eye, kron, block_diag
from scipy.sparse.linalg import splu
from scipy.linalg import lu_factor, lu_solve
from .regular_current import RegularCurrent
from .refined_kernel_induction import RefinedKernelInduction
from .variable_compatible_current import weak_vector
from .regular_fluid_linear import operator, periodic_derivative


class PlasmaPreconditioner:
    def __init__(self, coordinates):
        self.coordinates = c = coordinates; d = c.fluid; geo = c.geometry
        background = c.background(); model = RegularCurrent(c.system.mapping, background[..., 0], c.system.ion_mass)
        ids = c.fluid_indices; b_ids = c.magnetic_indices
        raw = operator(d, background)[np.ix_(ids, ids)]
        self.fluid = (diags(1/c.fluid_scale)@raw@diags(c.fluid_scale)).tocsc()
        solver = RefinedKernelInduction(model, **c.system.solver_options)
        currents, hall, lorentz = [], [], []
        bias = np.array([0., 0., c.bias])
        for k, index in enumerate(b_ids):
            magnetic = np.zeros(geo.hessian.shape[0]); magnetic[index] = c.magnetic_scale[k]
            current = model.current(magnetic); currents.append(current)
            force = np.cross(geo.sample_electric(current), bias)
            forcing = weak_vector(geo, model.delta*force/(E*model.density[..., None]))
            rate = solver.solve(forcing)
            hall.append(rate.magnetic[b_ids]/c.magnetic_scale)
            fluid = c.system.mapping.force_rate(background, force)
            lorentz.append(fluid.ravel()[ids]/c.fluid_scale)
        self.magnetic = np.column_stack(hall)
        self.to_fluid = np.column_stack(lorentz)
        current_matrix = np.column_stack(currents)
        inertia = np.column_stack([model.inertia@value for value in current_matrix.T])
        metric = geo.hessian[np.ix_(b_ids, b_ids)].toarray()*c.magnetic_scale[:, None]*c.magnetic_scale[None, :]
        metric += current_matrix.T@inertia
        # Only transverse momenta receive J cross B0. Their exact physical
        # energy Hessian gives the reciprocal ideal-induction work adjoint.
        kinetic = np.zeros(d.shape)
        kinetic[..., 1:3] = (d.volume*d.mass*d.r**2/c.density)[..., None]
        kinetic = kinetic.ravel()[ids]*c.fluid_scale**2
        self.to_magnetic = np.linalg.solve(metric, -self.to_fluid.T*kinetic[None, :])
        # Approximate uniform translation of magnetic coefficients; this is
        # preconditioning only, never the actual nonlinear current transport.
        mesh = geo.mesh; z = periodic_derivative(mesh.nz, mesh.dz, d.coefficients)
        advection = block_diag((kron(eye(mesh.nr+1), z), kron(eye(mesh.nr), z))).tocsr()
        self.magnetic += c.axial_speed*advection[np.ix_(b_ids, b_ids)].toarray()
        self._dt = None

    def factor(self, midpoint_interval):
        h = float(midpoint_interval)
        if not np.isfinite(h) or h <= 0:
            raise ValueError('Positive finite implicit interval required')
        if self._dt != h:
            self._fluid_factor = splu(eye(self.fluid.shape[0], format='csc')-h*self.fluid)
            self._fluid_response = self._fluid_factor.solve(self.to_fluid)
            schur = np.eye(self.magnetic.shape[0])-h*self.magnetic-h*h*self.to_magnetic@self._fluid_response
            self._magnetic_factor = lu_factor(schur)
            self._dt = h

    def solve(self, value):
        if self._dt is None:
            raise ValueError('Factor the implicit interval first')
        n = self.coordinates.nfluid; h = self._dt
        fluid = self._fluid_factor.solve(value[:n])
        magnetic = lu_solve(self._magnetic_factor, value[n:]+h*self.to_magnetic@fluid)
        fluid += h*self._fluid_response@magnetic
        return np.r_[fluid, magnetic]
