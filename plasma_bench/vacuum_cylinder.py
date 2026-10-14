"""Free plasma flux coupled to a magnetic-only periodic cylindrical exterior."""
from copy import copy
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.linalg import circulant
from scipy.sparse import csr_matrix, diags, kron, coo_matrix
from scipy.sparse.linalg import splu
from scipy.special import kve
from .compatible_geometry import gram
from .compatible_polynomials import bases
from .strong_wall_geometry import StrongWallGeometry
from .harmonic_magnetic import HarmonicMagnetic


class VacuumCylinderGeometry(StrongWallGeometry):
    def __init__(self, base, *, source_radius):
        if not np.isfinite(source_radius) or source_radius <= base.mesh.radius:
            raise ValueError('Ideal current-sheet source must be outside the plasma')
        self.__dict__.update(copy(base.__dict__))
        self.source_radius = float(source_radius)
        m = self.mesh; nr, nz = m.nr, m.nz
        self.theta_e = (self.Nx[:, 1:].multiply(1/self.radial_points[:, None])@diags(m.rf[1:])).tocsr()
        # Reuse full-current physical transport with an identity nodal transform.
        self.wall_projection = np.zeros(nr-1)
        self.pairs = ((self.radial_e, self.Nz), (self.theta_e, self.Nz), (self.Nx, self.Ez))
        self.mass = tuple(kron(gram(a, self.wr), gram(b, self.wz), format='csc') for a,b in self.pairs)
        self.solves = tuple(splu(a).solve for a in self.mass)
        self.offsets = np.cumsum([0]+[a.shape[0] for a in self.mass])
        old = np.r_[np.arange(base.offsets[2]), np.arange(base.offsets[2]+nz, self.offsets[-1])]
        embedding = csr_matrix((np.ones(len(old)), (old, np.arange(len(old)))),
                               shape=(self.offsets[-1], base.offsets[-1]))
        boundary = coo_matrix((-m.radius*np.ones(nz),
            (nr*nz+np.arange(nz), self.offsets[1]+(nr-1)*nz+np.arange(nz))),
            shape=(base.hessian.shape[0], self.offsets[-1])).tocsr()
        self.faraday = (base.faraday@embedding.T+boundary).tocsr()
        # Exact radial vacuum elimination for the SAME polynomial axial trace.
        mass_z = gram(self.Nz, self.wz)
        stiffness_z = gram(self.dNz, self.wz)
        zm = np.fft.fft(mass_z[:, 0].toarray().ravel()).real
        zs = np.fft.fft(stiffness_z[:, 0].toarray().ravel()).real
        if np.any(zm <= 0) or np.any(zs[1:] <= 0):
            raise ValueError('Positive nonconstant axial vacuum modes required')
        self.vacuum_k = np.r_[0., np.sqrt(zs[1:]/zm[1:])]
        k = self.vacuum_k[1:]; x = k*m.radius
        symbol = np.r_[0., 2*np.pi/MU0*zm[1:]*k*kve(0, x)/(m.radius*kve(1, x))]
        self.vacuum_metric = circulant(np.fft.ifft(symbol).real)
        rows, cols = np.meshgrid(nr*nz+np.arange(nz), nr*nz+np.arange(nz), indexing='ij')
        self.vacuum_hessian = coo_matrix((self.vacuum_metric.ravel(), (rows.ravel(), cols.ravel())),
                                        shape=base.hessian.shape).tocsc()
        self.plasma_hessian = base.hessian
        self.hessian = (base.hessian+self.vacuum_hessian).tocsc()

    def gradient(self, potential):
        p = np.asarray(potential)
        if p.shape != (self.mesh.nr+1, self.mesh.nz) or not np.isfinite(p).all():
            raise ValueError('Finite full nodal potential required')
        return np.r_[(self.Dx@p).ravel(), np.zeros(self.mesh.nr*self.mesh.nz), np.asarray(p@self.Dz.T).ravel()]

    def difference_curl(self, electric):
        m = self.mesh; r,t,z = self.split_electric(electric)
        r=r.reshape(m.nr,m.nz); t=t.reshape(m.nr,m.nz); z=z.reshape(m.nr+1,m.nz)
        p=np.vstack((np.zeros((1,m.nz)), -m.rf[1:,None]*t))
        toroidal=(r-np.roll(r,-1,axis=1))/m.dz+np.diff(z,axis=0)/m.radial_measure[:,None]
        return self.project_curl(np.r_[p.ravel(),toroidal.ravel()])

    def exterior_samples(self, wall_flux, radii, axial_points):
        """Induced Br,Btheta,Bz in the vacuum, with the actual polynomial trace."""
        m=self.mesh; r=np.asarray(radii,float); z=np.asarray(axial_points,float)
        if r.ndim != 1 or np.any(r<m.radius) or not np.isfinite(r).all():
            raise ValueError('Finite exterior radii required')
        wall=np.asarray(wall_flux,float)
        if wall.shape != (m.nz,) or not np.isfinite(wall).all():
            raise ValueError('Complete finite boundary flux required')
        n,dn,_,_=bases(m.zf,z,self.degree,periodic=True)
        phases=np.exp(2j*np.pi*np.arange(m.nz)[:,None]*np.arange(m.nz)[None,:]/m.nz)
        shapes=n@phases; derivatives=dn@phases
        amplitudes=np.fft.fft(wall)/m.nz
        k=self.vacuum_k[1:]; x=r[:,None]*k; root=m.radius*k
        decay=np.exp(-(x-root[None,:]))
        radial=kve(1,x)/kve(1,root)[None,:]*decay/m.radius
        axial=-k[None,:]*kve(0,x)/kve(1,root)[None,:]*decay/m.radius
        br=-(radial*amplitudes[None,1:])@derivatives[:,1:].T
        bz=(axial*amplitudes[None,1:])@shapes[:,1:].T
        return np.stack((br.real,np.zeros_like(br.real),bz.real),axis=-1)


class VacuumMagnetic(HarmonicMagnetic):
    """Affine source plus induced field with genuinely free outer flux."""
    def __init__(self, mesh, psi_increment, g, axial_bias):
        psi=np.asarray(psi_increment,float); g=np.asarray(g,float)
        if psi.shape != (mesh.nr+1,mesh.nz+1) or g.shape != (mesh.nr,mesh.nz):
            raise ValueError('Complete magnetic arrays required')
        if not np.isfinite(psi).all() or not np.isfinite(g).all() or not np.isfinite(axial_bias):
            raise ValueError('Finite field required')
        if np.any(psi[0] != 0) or np.any(psi[:,0] != psi[:,-1]):
            raise ValueError('Zero axis flux and periodic axial trace required')
        self.mesh=mesh; self.psi_increment=psi.copy(); self.g=g.copy(); self.axial_bias=float(axial_bias)
        self.psi_increment.setflags(write=False); self.g.setflags(write=False)

    def energy_components(self, model):
        geo=model.geometry; value=self.packed_increment(geo); length=2*self.mesh.half_length
        return dict(plasma_increment=float(.5*value@(geo.plasma_hessian@value)),
            vacuum_increment=float(.5*value@(geo.vacuum_hessian@value)),
            global_cross=float(2*np.pi*length*self.axial_bias*np.mean(self.psi_increment[-1,:-1])/MU0),
            source_background=float(np.pi*geo.source_radius**2*length*self.axial_bias**2/(2*MU0)))
