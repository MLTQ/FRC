"""Derivatives and one-sided physical traces of the frozen field definitions."""
import numpy as np
from scipy.sparse import coo_matrix,csr_matrix,diags
from .regular_field_map import RegularFieldMap,polynomial,periodic


def local_bases(edges,points,degree,periodic_axis=False,cells=None):
    count=len(edges)-1;n=count if periodic_axis else count+1
    rows=[];cols=[];val=[];der=[];er=[];ec=[];ev=[];ed=[]
    for row,point in enumerate(points):
        cell=min(np.searchsorted(edges,point,side='right')-1,count-1) if cells is None else cells[row]
        start=cell-(degree-1)//2
        if not periodic_axis:start=min(max(start,0),count-degree)
        ids=np.arange(start,start+degree+1);coords=edges[0]+ids*(edges[1]-edges[0]) if periodic_axis else edges[ids]
        scale=coords[-1]-coords[0];matrix=np.vander((coords-point)/scale,N=degree+1,increasing=True).T
        targets=np.zeros((degree+1,3));targets[0,0]=1;targets[1,1]=1/scale
        if degree>=2:targets[2,2]=2/scale**2
        value=np.linalg.solve(matrix,targets)
        rows.extend([row]*(degree+1));cols.extend(ids%n if periodic_axis else ids);val.extend(value[:,0]);der.extend(value[:,1])
        er.extend([row]*degree);ec.extend(ids[:-1]%count if periodic_axis else ids[:-1])
        ev.extend(-np.cumsum(value[:,1])[:-1]*np.diff(coords));ed.extend(-np.cumsum(value[:,2])[:-1]*np.diff(coords))
    sparse=lambda a,b,c,shape:coo_matrix((a,(b,c)),shape=shape).tocsr()
    return sparse(val,rows,cols,(len(points),n)),sparse(der,rows,cols,(len(points),n)),sparse(ev,er,ec,(len(points),count)),sparse(ed,er,ec,(len(points),count))


def current_maps(mesh,r,z,degree,radial_cells=None,axial_cells=None):
    nx,dx,ex,dex=local_bases(.5*mesh.rf**2,.5*r*r,degree,cells=radial_cells)
    nz,dz,ez,dez=local_bases(mesh.zf,z,degree,True,axial_cells)
    re=ex.multiply(r[:,None]);te=nx[:,1:-1].multiply(1/r[:,None])@diags(mesh.rf[1:-1])
    maps=((re,nz),(te,nz),(nx,ez))
    dr=((ex+dex.multiply(r[:,None]**2),nz),((dx[:,1:-1]-nx[:,1:-1].multiply(1/r[:,None]**2))@diags(mesh.rf[1:-1]),nz),(dx.multiply(r[:,None]),ez))
    dzmaps=((re,dz),(te,dz),(nx,dez))
    return maps,dr,dzmaps


class CurrentDifferential:
    def __init__(self,geometry,*,radial_cells=None,axial_cells=None,degree=None):
        self.geometry=geometry
        self.maps,self.dr,self.dz=current_maps(geometry.mesh,geometry.radial_points,geometry.axial_points,degree or geometry.degree,radial_cells,axial_cells)
        self.offsets=np.cumsum([0]+[a.shape[1]*b.shape[1] for a,b in self.maps])

    def sample(self,value,maps=None):
        maps=self.maps if maps is None else maps;g=self.geometry;result=[]
        for k,(a,b) in enumerate(maps):result.append(g.tensor(a,np.asarray(value)[self.offsets[k]:self.offsets[k+1]].reshape(a.shape[1],b.shape[1]),b))
        return np.stack(result,axis=-1)

    def adjoint(self,values,maps=None):
        maps=self.maps if maps is None else maps;g=self.geometry;result=[]
        for k,(a,b) in enumerate(maps):result.append(np.asarray((a.T@g._unlayout(g.weight*values[...,k]))@b).ravel())
        return np.concatenate(result)

    def gradient(self,value):
        v=self.sample(value);result=np.zeros(v.shape+(3,));result[...,0]=self.sample(value,self.dr);result[...,2]=self.sample(value,self.dz)
        result[...,0,1]=-v[...,1]/self.geometry.r;result[...,1,1]=v[...,0]/self.geometry.r
        return result

    def advective_adjoint(self,vector,velocity):
        result=self.adjoint(vector*velocity[...,0,None],self.dr)+self.adjoint(vector*velocity[...,2,None],self.dz)
        connection=np.zeros_like(vector);connection[...,0]=vector[...,1]*velocity[...,1]/self.geometry.r;connection[...,1]=-vector[...,0]*velocity[...,1]/self.geometry.r
        return result+self.adjoint(connection)


class TensorSurface:
    def __init__(self,geometry,radial):
        self.mesh=geometry.mesh;self.degree=geometry.degree
        self.radial_points=self.mesh.rf[1:-1] if radial else geometry.radial_points
        self.axial_points=geometry.axial_points if radial else self.mesh.zf[:-1]
        self.r,self.z=np.meshgrid(self.radial_points,self.axial_points,indexing='ij')
        wr=2*np.pi*self.radial_points if radial else geometry.wr;wz=geometry.wz if radial else np.ones(self.mesh.nz)
        self.weight=wr[:,None]*wz[None,:]

    def _unlayout(self,v):return np.asarray(v)
    def tensor(self,a,v,b):return np.asarray((a@v)@b.T)


def fluid_matrices(d,r,side=None,derivative=False,transverse=False):
    ratio=r*d.elements/d.radius;nearest=np.rint(ratio);at_face=abs(ratio-nearest)<32*np.finfo(float).eps*np.maximum(1,abs(ratio))
    element=np.floor(ratio).astype(int)
    if side is not None:element=np.where(at_face,nearest.astype(int)-(side=='left'),element)
    element=np.clip(element,0,d.elements-1);matrix=np.zeros((len(r),d.elements*(d.degree+1)))
    for e in range(d.elements):
        rows=np.flatnonzero(element==e);begin=1 if transverse and e==0 else 0;nodes=d.x[e,begin:,0];points=.5*r[rows]**2
        if derivative:
            values=np.zeros((len(rows),len(nodes)))
            for i in range(len(nodes)):
                for j in range(len(nodes)):
                    if i==j:continue
                    term=np.ones(len(rows))/(nodes[i]-nodes[j])
                    for k in range(len(nodes)):
                        if k!=i and k!=j:term*=(points-nodes[k])/(nodes[i]-nodes[k])
                    values[:,i]+=term
            values*=r[rows,None]
        else:values=polynomial(nodes,points)
        matrix[np.ix_(rows,e*(d.degree+1)+np.arange(begin,d.degree+1))]=values
    return csr_matrix(matrix)


def fluid_map(d,surface,side=None):
    mapping=RegularFieldMap(d,surface)
    mapping.radial=fluid_matrices(d,surface.radial_points,side);mapping.radial_velocity=fluid_matrices(d,surface.radial_points,side,transverse=True)
    return mapping


class RegularDifferential:
    def __init__(self,mapping):
        self.mapping=mapping;self.geometry=mapping.geometry;d=mapping.fluid
        self.dr=fluid_matrices(d,self.geometry.radial_points,derivative=True);self.drv=fluid_matrices(d,self.geometry.radial_points,derivative=True,transverse=True)
        theta=2*np.pi*(self.geometry.axial_points[:,None]-d.z.ravel()[None,:])/d.length;n=d.nz;matrix=np.zeros_like(theta)
        for k in range(1,(n-1)//2+1):matrix-=2*k*np.sin(k*theta)
        if n%2==0:matrix-=n/2*np.sin(n/2*theta)
        self.dz=csr_matrix(matrix*(2*np.pi/d.length/n))

    def scalar(self,value,axis,transverse=False):
        m=self.mapping;radial=(self.drv if transverse else self.dr) if axis==0 else (m.radial_velocity if transverse else m.radial)
        axial=self.dz if axis==2 else m.axial
        return self.geometry.tensor(radial,np.asarray(value).reshape(radial.shape[1],m.fluid.nz),axial)

    def adjoint(self,value,axis,transverse=False):
        m=self.mapping;radial=(self.drv if transverse else self.dr) if axis==0 else (m.radial_velocity if transverse else m.radial);axial=self.dz if axis==2 else m.axial
        return np.asarray((radial.T@self.geometry._unlayout(self.geometry.weight*value))@axial).reshape(m.shape)

    def velocity_gradient(self,state):
        m=self.mapping;g=self.geometry;v=m.velocity(state);n=m.scalar(state[...,0]);result=np.zeros(v.shape+(3,))
        for axis in (0,2):
            for k in (0,1):result[...,k,axis]=g.r*self.scalar(state[...,k+1]/state[...,0],axis,True)+(v[...,k]/g.r if axis==0 else 0)
            result[...,2,axis]=(self.scalar(state[...,3],axis)-v[...,2]*self.scalar(state[...,0],axis))/n
        result[...,0,1]=-v[...,1]/g.r;result[...,1,1]=v[...,0]/g.r
        return result

    def stress_rate(self,state,stress):
        m=self.mapping;g=self.geometry;d=m.fluid;n=m.scalar(state[...,0]);rate=np.zeros(d.shape)
        weak=[m.adjoint(stress[...,0,0]+stress[...,1,1],transverse=True)+self.adjoint(g.r*stress[...,0,0],0,True)+self.adjoint(g.r*stress[...,0,2],2,True),
              m.adjoint(stress[...,1,0]-stress[...,0,1],transverse=True)+self.adjoint(g.r*stress[...,1,0],0,True)+self.adjoint(g.r*stress[...,1,2],2,True)]
        den=d.mass*d.volume*2*d.x
        for k in (0,1):np.divide(weak[k],den,out=rate[...,k+1],where=den>0)
        wz=sum(self.adjoint(stress[...,2,axis]/n,axis)-m.adjoint(stress[...,2,axis]*self.scalar(state[...,0],axis)/n**2) for axis in (0,2))
        rate[...,3]=state[...,0]*wz/(d.mass*d.volume)
        return rate
