"""Original Gauss solver with numerical diffusion in its proposal inverse."""
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from scipy.linalg import block_diag
from .guarded_modal_gauss import GuardedModalGauss


class HyperdiffusiveModalGauss(GuardedModalGauss):
    def __init__(self, system, state, field, **options):
        super().__init__(system, state, field, **options)
        blocks = self.inverse.blocks
        c = self.coordinates
        mesh = c.geometry.mesh
        nr = mesh.nr
        for mode, operators in enumerate(self.inverse.operators):
            zn, ze, zd = (a[mode] for a in
                (blocks.node_symbol, blocks.edge_symbol, blocks.derivative_symbol))
            mass = block_diag(blocks.radial_e_mass*zn,
                blocks.theta_e_mass*zn, blocks.axial_e_mass*ze)
            hessian = block_diag((blocks.psi_radial*zd+blocks.psi_axial*zn)/MU0,
                blocks.radial_e_mass*ze/MU0)
            curl = np.zeros((2*nr-1, 3*nr), complex)
            curl[:nr-1, nr:2*nr-1] = -np.diag(mesh.rf[1:-1])
            curl[nr-1:, :nr] = np.eye(nr)*(1-np.exp(2j*np.pi*mode/c.nz))/mesh.dz
            curl[nr-1:, 2*nr-1:] = np.diff(np.eye(nr+1), axis=0)/mesh.radial_measure[:, None]
            current = -np.linalg.solve(mass, curl.conj().T@hessian@blocks.lift)
            h = blocks.lift.conj().T@hessian@blocks.lift
            d = current.conj().T@mass@current
            metric = h+blocks.reduced_mass/(E**2*c.density)*d
            laplacian = MU0*np.linalg.solve(h, d)
            damping = -system.eta6*np.linalg.solve(metric, d@laplacian@laplacian)
            fluid, to_fluid, to_magnetic, magnetic = operators
            self.inverse.operators[mode] = (fluid, to_fluid, to_magnetic, magnetic+damping)
