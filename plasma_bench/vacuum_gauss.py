"""Existing unsplit driven Gauss solve on the free-flux plasma/vacuum space."""
import numpy as np
from scipy.constants import elementary_charge as E
from scipy.linalg import block_diag
from scipy.sparse import block_diag as sparse_block_diag
from .modal_wall_coordinates import ModalWallCoordinates
from .modal_wall_blocks import ModalWallBlocks
from .modal_wall_inverse import ModalWallInverse
from .driven_modal_gauss import DrivenModalGauss
from .vacuum_cylinder import VacuumMagnetic
from .spitzer_wall_plasma import spitzer_samples


class VacuumCoordinates(ModalWallCoordinates):
    def __init__(self, system, state, field):
        super().__init__(system,state,field)
        self.magnetic_row_count += 1
        self.size += self.nz

    def pack(self, fluid, magnetic):
        old=super().pack(fluid,magnetic)
        psi,_=self.geometry.unpack(magnetic)
        insert=self.nfluid+(self.geometry.mesh.nr-1)*self.nz
        return np.r_[old[:insert],psi[-1,:-1]/self.magnetic_scalars[0],old[insert:]]

    def unpack(self, value):
        value=np.asarray(value,float)
        if value.shape != (self.size,) or not np.isfinite(value).all():
            raise ValueError('Finite complete vacuum/plasma coordinates required')
        insert=self.nfluid+(self.geometry.mesh.nr-1)*self.nz
        # Parent uses the declared size; retain its original contract through a
        # shallow coordinate view rather than mutate the live stage object.
        from copy import copy
        old=copy(self); old.size-=self.nz; old.magnetic_row_count-=1
        fluid,magnetic=ModalWallCoordinates.unpack(old,np.r_[value[:insert],value[insert+self.nz:]])
        psi,g=self.geometry.unpack(magnetic)
        psi[-1,:-1]=value[insert:insert+self.nz]*self.magnetic_scalars[0]; psi[-1,-1]=psi[-1,0]
        return fluid,self.geometry.pack(psi,g)


def symbols(matrix, nz):
    nr,nc=np.array(matrix.shape)//nz
    first=matrix[:,np.arange(nc)*nz].toarray().reshape(nr,nz,nc)
    return np.fft.fft(first,axis=1).transpose(1,0,2)


class VacuumBlocks(ModalWallBlocks):
    def __init__(self, coordinates, eta):
        super().__init__(coordinates)
        c=coordinates; g=c.geometry; nr=g.mesh.nr
        self.eta=eta
        self.mass_modes=symbols(sparse_block_diag(g.mass,format='csc'),c.nz)
        self.field_modes=symbols(g.hessian,c.nz); self.curl_modes=symbols(g.faraday,c.nz)
        self.lift=block_diag(np.vstack((np.zeros((1,nr)),np.eye(nr)))*c.magnetic_scalars[0],
                             c.radial_basis*c.magnetic_scalars[1])

    def coupled(self, mode):
        c=self.coordinates; g=c.geometry; nr=g.mesh.nr
        mass=self.mass_modes[mode]; hessian=self.field_modes[mode]; curl=self.curl_modes[mode]
        current=-np.linalg.solve(mass,curl.conj().T@hessian@self.lift)
        kinetic=current.conj().T@mass@current
        metric=self.lift.conj().T@hessian@self.lift+self.reduced_mass/(E**2*c.density)*kinetic
        hall=np.zeros_like(mass); hall[:nr,nr:2*nr]=self.cross_e*self.node_symbol[mode]
        hall[nr:2*nr,:nr]=-self.cross_e.T*self.node_symbol[mode]
        magnetic=np.linalg.solve(metric,-self.eta*kinetic-self.delta*c.bias/(E*c.density)*current.conj().T@hall@current)
        magnetic+=np.eye(magnetic.shape[0])*c.axial_speed*self.z_symbol[mode]
        phase=np.exp(2j*np.pi*mode*np.arange(c.nz)/c.nz)
        cross_z=np.sum(g.wz*(c.system.mapping.axial@phase).conj()*(g.Nz@phase))/c.nz
        force=np.zeros((len(self.force_denominator),6,3*nr+1),complex); active=self.force_denominator>0
        force[active,1,nr:2*nr]=c.bias*cross_z*self.force_theta[active]/self.force_denominator[active,None]
        force[active,2,:nr]=-c.bias*cross_z*self.force_radial[active]/self.force_denominator[active,None]
        force=force.reshape(-1,3*nr+1)[c.fluid_rows]/c.fluid_row_scale[:,None]
        to_fluid=force@current
        to_magnetic=np.linalg.solve(metric,-to_fluid.conj().T*self.kinetic[None,:])
        return to_fluid,to_magnetic,magnetic


class VacuumInverse(ModalWallInverse):
    def __init__(self, coordinates, eta):
        self.coordinates=coordinates; self.blocks=VacuumBlocks(coordinates,eta)
        self.operators=[(self.blocks.fluid(k),*self.blocks.coupled(k)) for k in range(coordinates.nz)]
        self._shift=None


class VacuumGauss(DrivenModalGauss):
    field_type=VacuumMagnetic

    def __init__(self, system, state, field, *, absolute_tolerance=1e-12,
                 relative_tolerance=1e-10, max_iterations=60):
        if not np.isfinite([absolute_tolerance,relative_tolerance]).all() or min(absolute_tolerance,relative_tolerance)<=0:
            raise ValueError('Positive finite nonlinear tolerances required')
        if type(max_iterations) is not int or max_iterations<1:raise ValueError('Positive iteration budget required')
        self.system=system; self.coordinates=VacuumCoordinates(system,state,field)
        eta=system.resistivity if system.transport=='fixed' else float(np.mean(spitzer_samples(system.mapping,state,system.ion_mass)[2]))
        self.inverse=VacuumInverse(self.coordinates,eta)
        self.atol=absolute_tolerance; self.rtol=relative_tolerance; self.max_iterations=max_iterations
        root=np.sqrt(3.)/6; self.a=np.array([[.25,.25-root],[.25+root,.25]])
