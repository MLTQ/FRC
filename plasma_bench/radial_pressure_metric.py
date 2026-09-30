"""Diagnostic radial pressure force paired with the true variable-density kinetic metric."""
import numpy as np
from scipy.sparse import kron
from scipy.sparse.linalg import LinearOperator,cg,splu
from .compatible_geometry import gram


class RadialPressureMetric:
    def __init__(self,kinematics,density_map):
        self.kinematics=kinematics;self.geometry=kinematics.geometry;self.mesh=kinematics.mesh
        self.radial=kinematics.radial[1];self.axial=density_map.axial
        g=self.geometry
        self.factor=splu(kron(gram(self.radial,g.wr),gram(self.axial,g.wz),format='csc'))

    def solve(self,model,flux,work):
        g=self.geometry;shape=(self.mesh.nr,self.mesh.nz);n=model.density
        scale=float(np.sum(g.weight*n)/self.mesh.volumes.sum());inverse=scale/n
        def action(flat):
            field=g.tensor(self.radial,np.asarray(flat).reshape(shape),self.axial)
            return np.asarray(self.radial.T@g._unlayout(g.weight*inverse*field)@self.axial).ravel()
        load=-scale/self.kinematics.mass*np.asarray(flux.radial_face.T@(self.mesh.rf[:,None]*work))
        if not np.any(load):return np.zeros(shape),0.
        operator=LinearOperator((load.size,load.size),matvec=action,dtype=float)
        preconditioner=LinearOperator(operator.shape,matvec=self.factor.solve,dtype=float)
        value,info=cg(operator,load.ravel(),M=preconditioner,rtol=1e-13,atol=0.,maxiter=500)
        residual=float(np.linalg.norm(action(value)-load.ravel())/np.linalg.norm(load))
        if info or not np.isfinite(value).all() or residual>2e-12:
            raise RuntimeError(f'Radial physical pressure metric solve failed: {info}, {residual}')
        return value.reshape(shape),residual
