"""Fine material density cuts with a separate continuous velocity space."""
import numpy as np
from scipy.constants import elementary_charge as EV
from scipy.linalg import eigh, solve_banded
from .species_material import MaterialSpecies, SpeciesMaterialState


class MaterialSubcells(MaterialSpecies):
    """One physical species; every density cut retained, nested P1 mechanics."""
    def __init__(self, numbers, mass, charge, bounds, faces, velocity_cells, *, number_scale=None, **controls):
        numbers=np.asarray(numbers,float)
        self.number_scale=float(numbers.max()) if number_scale is None else float(number_scale)
        if not np.isfinite(self.number_scale) or self.number_scale<=0:raise ValueError('Positive finite population coordinate scale required')
        self.relative_numbers=numbers/self.number_scale if number_scale is None else numbers.copy()
        super().__init__(self.relative_numbers[None,:], [mass], [charge], bounds, **controls)
        # Physical products are diagnostic only. All local evolution retains
        # the factored pair even when their product is below float64 range.
        self.numbers=(self.relative_numbers*self.number_scale)[None,:]
        if not isinstance(velocity_cells,int) or velocity_cells<1:raise ValueError('Positive velocity cell count required')
        self.reference_faces=np.asarray(faces,np.longdouble).copy()
        if self.reference_faces.shape!=(self.cells+1,) or np.any(np.diff(self.reference_faces)<=0):raise ValueError('Ordered physical density faces required')
        target=self.reference_faces[0]+np.ptp(self.reference_faces)*(1-np.cos(np.linspace(0,np.pi,velocity_cells+1)))/2
        self.velocity_indices=np.unique(np.argmin(abs(self.reference_faces[:,None]-target[None,:]),axis=0))
        self.velocity_cells=len(self.velocity_indices)-1
        coarse=self.reference_faces[self.velocity_indices]
        self.owner=np.searchsorted(coarse,self.reference_faces[:-1],side='right')-1
        width=np.diff(coarse)[self.owner]
        self.span=np.asarray(np.diff(self.reference_faces)/width,float)
        self.left_shape=np.asarray((coarse[self.owner+1,None]-np.stack((self.reference_faces[:-1],self.reference_faces[1:]),axis=1))/width[:,None],float)
        self.right_shape=np.asarray((np.stack((self.reference_faces[:-1],self.reference_faces[1:]),axis=1)-coarse[self.owner,None])/width[:,None],float)
        l,r=self.left_shape,self.right_shape;cell_mass=self.mass[0]*self.relative_numbers
        ll=(l[:,0]**2+l[:,0]*l[:,1]+l[:,1]**2)/3
        rr=(r[:,0]**2+r[:,0]*r[:,1]+r[:,1]**2)/3
        lr=(2*l[:,0]*r[:,0]+l[:,0]*r[:,1]+l[:,1]*r[:,0]+2*l[:,1]*r[:,1])/6
        matrix=np.zeros((self.velocity_cells+1,self.velocity_cells+1))
        np.add.at(matrix,(self.owner,self.owner),cell_mass*ll)
        np.add.at(matrix,(self.owner+1,self.owner+1),cell_mass*rr)
        np.add.at(matrix,(self.owner,self.owner+1),cell_mass*lr)
        np.add.at(matrix,(self.owner+1,self.owner),cell_mass*lr)
        self.mass_matrix=matrix;self.mass_scale=float(np.diag(matrix).max())
        self.matrix=np.zeros((1,3,self.velocity_cells+1));self.matrix[0,1]=np.diag(matrix)
        self.matrix[0,0,1:]=np.diag(matrix,1);self.matrix[0,2,:-1]=np.diag(matrix,-1)

    def validate(self,state):
        if (state.faces.shape!=(1,self.cells+1) or state.velocity.shape!=(1,self.velocity_cells+1,3) or
                state.thermal.shape!=(1,self.cells) or not all(np.isfinite(getattr(state,k)).all() for k in state.__dataclass_fields__)):
            raise ValueError('Complete subcell geometry, coarse velocity and fine thermal state required')
        width=np.diff(state.faces,axis=1)
        if np.any(width<=0) or state.faces[0,0]<=self.bounds[0] or state.faces[0,-1]>=self.bounds[1]:raise ValueError('Density subcells crossed or reached remote boundary')
        if np.any(state.thermal<0):raise ValueError('Negative subcell thermal energy; no repair')
        return width

    def state(self,faces,temperature_ev,velocity=None):
        velocity=np.zeros((1,self.velocity_cells+1,3)) if velocity is None else np.asarray(velocity,float).copy()
        result=SpeciesMaterialState(np.asarray(faces,np.longdouble).copy(),velocity,self.relative_numbers[None,:]*EV*np.broadcast_to(temperature_ev,self.numbers.shape)/(self.gamma-1))
        self.validate(result);return result

    def fine_velocity(self,velocity):
        v=velocity[0];j=self.owner
        left=self.left_shape[:,0,None]*v[j]+self.right_shape[:,0,None]*v[j+1]
        right=self.left_shape[:,1,None]*v[j]+self.right_shape[:,1,None]*v[j+1]
        return left,right

    def restrict_load(self,load):
        # Coarse shapes are exactly linear combinations of nested fine shapes.
        result=np.zeros((self.velocity_cells+1,3));j=self.owner
        np.add.at(result,j,self.left_shape[:,0,None]*load[:-1])
        np.add.at(result,j+1,self.right_shape[:,0,None]*load[:-1])
        result[-1]+=load[-1]
        return result

    def pressure(self,state):
        width=self.validate(state);pressure=(self.gamma-1)*state.thermal/width
        sound=np.sqrt(self.gamma*pressure*width/(self.mass[0]*self.relative_numbers[None,:]))
        coarse=state.faces[0,self.velocity_indices];h=np.diff(coarse);jump=np.diff(state.velocity[0,:,2]);compression=np.maximum(-jump,0.)
        if self.viscosity_limiter=='minmod':
            gradient=jump/h;left=np.r_[gradient[0],gradient[:-1]];right=np.r_[gradient[1:],gradient[-1]]
            resolved=np.minimum(np.minimum(np.maximum(-left,0.),np.maximum(-gradient,0.)),np.maximum(-right,0.))
            compression=np.maximum(compression-h*resolved,0.)
        c=compression[self.owner][None,:]
        viscosity=self.mass[0]*self.relative_numbers[None,:]/width*(self.viscosity[0]*sound*c+self.viscosity[1]*c*c)
        return pressure,viscosity,sound

    def physical_thermal(self,state):
        return state.thermal*self.number_scale

    def number_density(self,state):
        return (self.relative_numbers[None,:]/self.validate(state))*self.number_scale

    def acoustic_timestep(self,state):
        p,q,_=self.pressure(state);h=np.diff(state.faces[0,self.velocity_indices]);width=self.validate(state)[0]
        coefficient=np.asarray(self.gamma*(p[0]+q[0])*width/h[self.owner]**2,float)
        stiffness=np.zeros_like(self.mass_matrix);j=self.owner
        np.add.at(stiffness,(j,j),coefficient);np.add.at(stiffness,(j+1,j+1),coefficient)
        np.add.at(stiffness,(j,j+1),-coefficient);np.add.at(stiffness,(j+1,j),-coefficient)
        frequency=0.
        if np.any(coefficient):
            eigen=eigh(stiffness/self.mass_scale,self.mass_matrix/self.mass_scale,
                subset_by_index=[self.velocity_cells,self.velocity_cells],eigvals_only=True,check_finite=False)[0]
            if eigen<0:raise ValueError('Negative maximum acoustic eigenvalue')
            frequency=np.sqrt(eigen)
        jump=abs(np.diff(state.velocity[0,:,2]));deformation=float(np.max(jump/h,initial=0.))
        return self.cfl/max(frequency,deformation) if max(frequency,deformation)>0 else np.inf

    def material_rate(self,state,geometry,frame_z):
        p,q,_=self.pressure(state);stress=p+q;j=self.owner;force=np.zeros_like(state.velocity)
        weighted=np.asarray(stress[0]*self.span,float)
        np.add.at(force[0,:,2],j,-weighted);np.add.at(force[0,:,2],j+1,weighted)
        edges,electric,active,_=geometry;cell,inside=active[0];k=cell[inside];owner=self.owner[k]
        h=np.diff(edges)[inside];coarse=state.faces[0,self.velocity_indices];width=np.diff(state.faces[0])[k]
        nodes=np.array([.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)])
        right=((edges[:-1][inside]-coarse[owner])[:,None]+h[:,None]*nodes)/np.diff(coarse)[owner,None]
        field=electric[:-1][inside,None]+np.diff(electric)[inside,None]*nodes
        weight=EV*self.charge[0]*self.relative_numbers[k]/width*h/2
        np.add.at(force[0,:,2],owner,weight*np.sum(field*(1-right),axis=1))
        np.add.at(force[0,:,2],owner+1,weight*np.sum(field*right,axis=1))
        acceleration=solve_banded((1,1),self.matrix[0]/self.mass_scale,force[0]/self.mass_scale,check_finite=False)[None,...]
        left,right_velocity=self.fine_velocity(state.velocity)
        face_velocity=np.r_[left[:,2],right_velocity[-1,2]][None,:]+frame_z
        jump=self.span*np.diff(state.velocity[0,:,2])[j]
        return SpeciesMaterialState(face_velocity,acceleration,-stress*jump[None,:])

    def sample(self,state,points):
        x=np.asarray(points);width=self.validate(state)[0];cell=np.searchsorted(state.faces[0],x,side='right')-1
        inside=(cell>=0)&(cell<self.cells);j=cell[inside];owner=self.owner[j];coarse=state.faces[0,self.velocity_indices]
        a=(x[inside]-coarse[owner])/np.diff(coarse)[owner]
        number=np.zeros((1,len(x)));velocity=np.zeros((1,len(x),3));temperature=np.zeros_like(number)
        number[0,inside]=(self.relative_numbers[j]/width[j])*self.number_scale
        velocity[0,inside]=(1-a[:,None])*state.velocity[0,owner]+a[:,None]*state.velocity[0,owner+1]
        temperature[0,inside]=(self.gamma-1)*state.thermal[0,j]/(self.relative_numbers[j]*EV)
        # MaterialSpeciesSet supplies the common field after all blocks gather.
        return dict(number=number,velocity=velocity,temperature_ev=temperature)
