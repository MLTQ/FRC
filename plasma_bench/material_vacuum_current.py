"""Couple the full moving finite-mass plasma to its magnetic-only exterior."""
import numpy as np
from scipy.constants import elementary_charge as E,mu_0 as MU0
from scipy.sparse import csr_matrix
from .material_vector_current import MaterialVectorCurrent
from .material_vector_geometry import vector_geometry
from .material_carrier_current import MaterialCurrentState,positive_solve


class MaterialVacuumCurrent(MaterialVectorCurrent):
    def __init__(self,fluid,space,vacuum,**options):
        super().__init__(fluid,space,**options);self.vacuum=vacuum
        if vacuum.degree!=space.degree or space.vacuum_trace.shape[0]!=vacuum.boundary_count:
            raise ValueError('Identical physical polynomial boundary trace required')

    def geometry(self,state,fields=None):
        if fields is not None:raise ValueError('One complete moving geometry required')
        g=vector_geometry(self,state);v=self.vacuum.geometry(state.coordinates,g['fields']['velocity'])
        trace=self.space.vacuum_trace.toarray();g['h']+=trace.T@v['hessian']@trace
        g['hdot']+=trace.T@v['hessian_rate']@trace
        nt=self.space.toroidal_magnetic_count;ne=self.space.meridional_electric_count
        # A theta electric coefficient curls into the corresponding complete
        # poloidal magnetic coefficient, including its harmonic vacuum extension.
        d=g['d'].toarray();ddot=g['ddot'].toarray();d[:,ne:]=-g['h'][:,nt:];ddot[:,ne:]=-g['hdot'][:,nt:]
        g.update(d=csr_matrix(d),ddot=csr_matrix(ddot),vacuum_geometry=v)
        return g

    def project_initial(self,fluid_state,magnetic_field,current_density,boundary_flux):
        g=self.geometry(fluid_state);b=np.asarray(magnetic_field,float).reshape(-1,3);j=np.asarray(current_density,float).reshape(-1,3)
        boundary=np.asarray(boundary_flux,float)
        if b.shape!=(len(self.particles),3) or j.shape!=b.shape or boundary.shape!=(self.vacuum.boundary_count,):raise ValueError('Consistent physical interior field/current and exterior trace required')
        load=sum(g['magnetic'][i].T@(g['volume']/MU0*b[:,i]) for i in range(3))
        load+=self.space.vacuum_trace.T@(g['vacuum_geometry']['hessian']@boundary)
        drift=sum(g['electric'][i].T@(g['volume']/E*j[:,i]) for i in range(3))
        operator=g['h']+self.mu/E**2*(g['d']@g['solve'](g['d'].T.toarray()))
        rhs=load-self.mu/E*(g['d']@g['solve'](drift))
        return MaterialCurrentState(fluid_state,positive_solve(operator,np.asarray(rhs)))

    def response(self,state,*,rate=True):
        a=super().response(state,rate=rate)
        vacuum=self.vacuum.field_response(a['vacuum_geometry'],self.space.vacuum_trace@state.magnetic)
        a.update(vacuum_response=vacuum,vacuum_magnetic_energy=vacuum['energy'])
        if rate:
            work=float(np.sum(vacuum['force']*a['fields']['velocity']))
            a['force']+=vacuum['force'];a['mechanical_power']+=work;a['pair_power_residual']+=work
        return a
