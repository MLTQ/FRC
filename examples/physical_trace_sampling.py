"""Exact one-sided polynomial traces of actual fluid and electromagnetic fields."""
import numpy as np
from numpy.polynomial import Polynomial
from .current_transport_sampling import TransportSampler


def selected_bases(edges,points,cells,degree,periodic=False):
    count=len(edges)-1;size=count if periodic else count+1
    nodal=[np.zeros((len(points),size)) for _ in range(3)]
    edge=[np.zeros((len(points),count)) for _ in range(2)]
    for row,(point,cell) in enumerate(zip(points,cells)):
        start=int(cell)-(degree-1)//2
        if not periodic:start=min(max(start,0),count-degree)
        ids=np.arange(start,start+degree+1)
        coords=edges[0]+ids*(edges[1]-edges[0]) if periodic else edges[ids]
        scale=coords[-1]-coords[0];roots=(coords-point)/scale;derivatives=[]
        for i,index in enumerate(ids):
            others=np.delete(roots,i);poly=Polynomial.fromroots(others)/np.prod(roots[i]-others)
            values=[poly.deriv(k)(0.)/scale**k for k in range(3)]
            for k in range(3):nodal[k][row,index%size]+=values[k]
            derivatives.append(values)
        for i,index in enumerate(ids[:-1]):
            for k in range(2):edge[k][row,index%count]=-sum(v[k+1] for v in derivatives[:i+1])*(coords[i+1]-coords[i])
    return *nodal,*edge


class TraceSampler(TransportSampler):
    def __init__(self,fluid,mesh,r,z,radial_cells,axial_cells):
        self.d=fluid;self.mesh=mesh;self.r=np.asarray(r);self.z=np.asarray(z)
        radial_cells=np.asarray(radial_cells);axial_cells=np.asarray(axial_cells)
        if mesh.nr%fluid.elements or self.r.shape!=radial_cells.shape or self.z.shape!=axial_cells.shape:
            raise ValueError('Aligned meshes and a selected cell for every point required')
        if np.any(self.r<=0) or np.any(radial_cells<0) or np.any(radial_cells>=mesh.nr):
            raise ValueError('Positive radii and interior radial cells required; axis face has zero area')
        self.R=self.r[:,None];self.Z=self.z[None,:]
        self.nx,self.dx,self.dxx,self.ex,self.dex=selected_bases(.5*mesh.rf**2,.5*self.r**2,radial_cells,5)
        self.nz,self.dz,self.dzz,self.ez,self.dez=selected_bases(mesh.zf,self.z,axial_cells,5,True)
        d=fluid;phase=(self.z+.5*d.length)/d.length;freq=np.fft.fftfreq(d.nz)*d.nz
        self.fourier=np.exp(2j*np.pi*phase[:,None]*freq[None,:]);self.fourier_d=self.fourier*(2j*np.pi*freq/d.length)
        if d.nz%2==0:
            self.fourier[:,d.nz//2]=np.cos(np.pi*d.nz*phase)
            self.fourier_d[:,d.nz//2]=-np.pi*d.nz/d.length*np.sin(np.pi*d.nz*phase)
        self.radial={};elements=radial_cells*d.elements//mesh.nr
        for transverse in (False,True):
            value=np.zeros((len(self.r),d.elements*(d.degree+1)));derivative=value.copy()
            for e in range(d.elements):
                rows=np.flatnonzero(elements==e);begin=1 if transverse and e==0 else 0
                nodes=d.x[e,begin:,0];scale=nodes[-1]-nodes[0];origin=nodes[0]
                roots=(nodes-origin)/scale;points=(.5*self.r[rows]**2-origin)/scale
                for i in range(len(nodes)):
                    others=np.delete(roots,i);factors=points[:,None]-others[None,:]
                    denominator=np.prod(roots[i]-others);col=e*(d.degree+1)+begin+i
                    value[rows,col]=np.prod(factors,axis=1)/denominator
                    derivative[rows,col]=sum(np.prod(np.delete(factors,j,axis=1),axis=1)
                        for j in range(len(others)))/denominator*self.r[rows]/scale
            self.radial[transverse]=(value,derivative)
