"""Azimuthal current space including its independent driven-wall trace."""
from copy import copy
import numpy as np
from scipy.sparse import hstack, diags, kron, eye, csr_matrix
from scipy.sparse.linalg import splu
from .compatible_geometry import gram
from .strong_wall_geometry import StrongWallGeometry
from .physical_derivatives import CurrentDifferential, RegularDifferential, TensorSurface, fluid_map, local_bases
from .current_transport import CurrentTransport


class WallCurrentGeometry(StrongWallGeometry):
    def __init__(self, base):
        self.__dict__.update(copy(base.__dict__))
        self.base_geometry = base
        m = self.mesh; r = self.radial_points
        wall = base.Nx[:, -1].multiply(m.radius/r[:, None]).tocsr()
        radial_mass = gram(base.theta_e, self.wr)
        cross = (base.theta_e.T@wall.multiply(self.wr[:, None])).toarray().ravel()
        self.wall_projection = splu(radial_mass).solve(cross)
        # chi(R)=1; chi is orthogonal to the original interior theta basis.
        chi = wall-csr_matrix((base.theta_e@self.wall_projection)[:, None])
        self.theta_e = hstack((base.theta_e, chi), format='csr')
        self.pairs = ((self.radial_e, self.Nz), (self.theta_e, self.Nz), (self.Nx, self.Ez))
        self.mass = tuple(kron(gram(a, self.wr), gram(b, self.wz), format='csc') for a,b in self.pairs)
        self.solves = tuple(splu(a).solve for a in self.mass)
        self.offsets = np.cumsum([0]+[a.shape[0] for a in self.mass])
        begin = base.offsets[2]
        self.wall_indices = np.arange(begin, begin+m.nz)
        self.old_indices = np.r_[np.arange(begin), np.arange(begin+m.nz, self.offsets[-1])]
        self.embedding = csr_matrix((np.ones(len(self.old_indices)),
            (self.old_indices, np.arange(len(self.old_indices)))), shape=(self.offsets[-1],base.offsets[-1]))
        self.faraday = (base.faraday@self.embedding.T).tocsr()

    def embed(self, old, wall=None):
        value = np.zeros(self.offsets[-1]); value[self.old_indices] = old
        if wall is not None:
            if np.shape(wall) != (self.mesh.nz,) or not np.isfinite(wall).all():
                raise ValueError('Finite complete wall-current vector required')
            value[self.wall_indices] = wall
        return value

    def gradient(self, potential):
        return self.embed(self.base_geometry.gradient(potential))

    def nodal_theta(self, current):
        """Physical azimuthal nodal values, including the axis and wall."""
        nr,nz=self.mesh.nr,self.mesh.nz
        theta=self.split_electric(current)[1].reshape(nr,nz)
        return np.vstack((np.zeros((1,nz)),theta[:-1]-self.wall_projection[:,None]*theta[-1],theta[-1:]))


class WallCurrentDifferential(CurrentDifferential):
    def __init__(self, geometry, projection, *, radial_cells=None, axial_cells=None):
        self.geometry=geometry; m=geometry.mesh; r=geometry.radial_points
        nx,dx,ex,dex=local_bases(.5*m.rf**2,.5*r*r,geometry.degree,cells=radial_cells)
        nz,dz,ez,dez=local_bases(m.zf,geometry.axial_points,geometry.degree,True,axial_cells)
        transform=np.eye(m.nr);transform[:-1,-1]=-projection
        nodal=(nx[:,1:]@diags(m.rf[1:])@csr_matrix(transform)).tocsr()
        deriv=(dx[:,1:]@diags(m.rf[1:])@csr_matrix(transform)).tocsr()
        re=ex.multiply(r[:,None]);te=nodal.multiply(1/r[:,None])
        self.maps=((re,nz),(te,nz),(nx,ez))
        self.dr=((ex+dex.multiply(r[:,None]**2),nz),
            (deriv-nodal.multiply(1/r[:,None]**2),nz),(dx.multiply(r[:,None]),ez))
        self.dz=((re,dz),(te,dz),(nx,dez))
        self.offsets=np.cumsum([0]+[a.shape[1]*b.shape[1] for a,b in self.maps])


class WallCurrentTransport(CurrentTransport):
    """Original full transport pair on the enriched physical current fields."""
    def __init__(self, model):
        self.model=model;g=model.geometry;p=g.wall_projection
        self.current=WallCurrentDifferential(g,p);self.fluid=RegularDifferential(model.mapping);self.faces=[]
        for radial in (True,False):
            surface=TensorSurface(g,radial);count=g.mesh.nr if radial else g.mesh.nz
            right=np.arange(1,count) if radial else np.arange(count);left=right-1
            choices=({'radial_cells':left},{'radial_cells':right}) if radial else ({'axial_cells':left},{'axial_cells':right})
            currents=[WallCurrentDifferential(surface,p,**choice) for choice in choices]
            fluids=[fluid_map(model.mapping.fluid,surface,side if radial else None) for side in ('left','right')]
            self.faces.append((radial,surface,currents,fluids))
