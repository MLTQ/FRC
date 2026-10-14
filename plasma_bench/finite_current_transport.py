"""Physical current transport on a finite side/end field space."""
import numpy as np
from scipy.sparse import diags
from .physical_derivatives import CurrentDifferential,TensorSurface,local_bases
from .current_transport import CurrentTransport
from .finite_cylinder_map import FiniteCylinderMap,FiniteFluidDifferential


class FiniteCurrentDifferential(CurrentDifferential):
    def __init__(self,geometry,*,radial_cells=None,axial_cells=None):
        self.geometry=geometry;m=geometry.mesh;r=geometry.radial_points
        nx,dx,ex,dex=local_bases(.5*m.rf**2,.5*r*r,geometry.degree,cells=radial_cells)
        nz,dz,ez,dez=local_bases(m.zf,geometry.axial_points,geometry.degree,cells=axial_cells)
        nodal=nx[:,1:]@diags(m.rf[1:]);derivative=dx[:,1:]@diags(m.rf[1:])
        re=ex.multiply(r[:,None]);te=nodal.multiply(1/r[:,None])
        self.maps=((re,nz),(te,nz),(nx,ez))
        self.dr=((ex+dex.multiply(r[:,None]**2),nz),(derivative-nodal.multiply(1/r[:,None]**2),nz),(dx.multiply(r[:,None]),ez))
        self.dz=((re,dz),(te,dz),(nx,dez))
        self.offsets=np.cumsum([0]+[a.shape[1]*b.shape[1] for a,b in self.maps])


class FiniteSurface(TensorSurface):
    def __init__(self,geometry,radial):
        self.mesh=geometry.mesh;self.degree=geometry.degree
        self.radial_points=self.mesh.rf[1:-1] if radial else geometry.radial_points
        self.axial_points=geometry.axial_points if radial else self.mesh.zf[1:-1]
        self.r,self.z=np.meshgrid(self.radial_points,self.axial_points,indexing='ij')
        wr=2*np.pi*self.radial_points if radial else geometry.wr
        wz=geometry.wz if radial else np.ones(self.mesh.nz-1)
        self.weight=wr[:,None]*wz[None,:]


class FiniteCurrentTransport(CurrentTransport):
    def __init__(self,model):
        self.model=model;g=model.geometry
        self.current=FiniteCurrentDifferential(g);self.fluid=FiniteFluidDifferential(model.mapping);self.faces=[]
        for radial in (True,False):
            surface=FiniteSurface(g,radial);count=g.mesh.nr if radial else g.mesh.nz
            right=np.arange(1,count);left=right-1
            choices=({'radial_cells':left},{'radial_cells':right}) if radial else ({'axial_cells':left},{'axial_cells':right})
            currents=[FiniteCurrentDifferential(surface,**choice) for choice in choices]
            fluids=[FiniteCylinderMap(model.mapping.fluid,surface,side=side if radial else None) for side in ('left','right')]
            self.faces.append((radial,surface,currents,fluids))
