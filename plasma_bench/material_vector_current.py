"""Moving full-vector finite-mass plasma with species transport and heat."""
import numpy as np
from scipy.constants import elementary_charge as E,mu_0 as MU0
from .material_carrier_current import MaterialCarrierCurrent,MaterialCurrentState,positive_solve
from .material_vector_geometry import vector_geometry
from .material_dual_entropy import MaterialDualEntropy
from .material_collision_sources import collision_config,collision_sources


class MaterialVectorCurrent(MaterialCarrierCurrent):
    """All axisymmetric field/current/flow components; compact magnetic boundary."""
    def __init__(self,fluid,space,*,resistivity=.001,coulomb_exchange=True,collision_model=None):
        super().__init__(fluid,space)
        if not np.isfinite(resistivity) or resistivity<0:raise ValueError('Finite nonnegative prescribed resistivity required')
        self.resistivity=float(resistivity);self.coulomb_exchange=bool(coulomb_exchange)
        self.collision_model=collision_config(collision_model,fluid,self.resistivity,self.coulomb_exchange)
        self.electron_heat_fraction=fluid.ion_mass/fluid.mass
        self.species=MaterialDualEntropy(fluid,space)

    def geometry(self,state,fields=None):
        if fields is not None:raise ValueError('Full vector geometry evaluates one actual complete fluid state')
        return vector_geometry(self,state)

    def project_initial(self,fluid_state,magnetic_field,current_density):
        g=self.geometry(fluid_state);b=np.asarray(magnetic_field,float).reshape(-1,3);j=np.asarray(current_density,float).reshape(-1,3)
        if b.shape!=(len(self.particles),3) or j.shape!=b.shape or not np.isfinite(b).all() or not np.isfinite(j).all():
            raise ValueError('Consistent finite analytic vector B and J at material quadrature required')
        magnetic_load=sum(g['magnetic'][i].T@(g['volume']/MU0*b[:,i]) for i in range(3))
        drift_load=sum(g['electric'][i].T@(g['volume']/E*j[:,i]) for i in range(3))
        operator=g['h']+self.mu/E**2*(g['d']@g['solve'](g['d'].T.toarray()))
        rhs=magnetic_load-self.mu/E*(g['d']@g['solve'](drift_load))
        return MaterialCurrentState(fluid_state,positive_solve(operator,np.asarray(rhs)))

    def response(self,state,*,rate=True):
        f=self.fluid;g=self.geometry(state.fluid);b=state.magnetic
        if b.shape!=(self.space.magnetic_count,) or not np.isfinite(b).all():raise ValueError('Finite complete vector magnetic state required')
        h=g['h'];d=g['d'];mass=g['mass'];solve=g['solve'];w=solve(-d.T@b/E)
        drift=np.column_stack([v@w for v in g['electric']]);magnetic=np.column_stack([v@b for v in g['magnetic']])
        wb=float(.5*b@h@b);kj=float(.5*self.mu*w@(mass@w))
        result=dict(g,drift_coefficients=w,drift=drift,magnetic_field=magnetic,magnetic_energy=wb,current_energy=kj,
            total_energy=g['fields']['energy']+wb+kj,poloidal_flux=self.space.flux@b)
        if not rate:return result
        thermal=self.species.pair(state.fluid,result,drift)
        curlw=np.column_stack([c@w for c in g['curl']])
        physical_forcing=self.delta*np.cross(drift,magnetic+self.mu/E*curlw)-self.mu/E*np.cross(drift,g['vorticity'])
        n=g['fields']['n'].ravel()
        collisions=collision_sources(g['fields'],drift,f.ion_mass,model=self.collision_model,
            resistivity=self.resistivity,coulomb_exchange=self.coulomb_exchange)
        physical_forcing+=collisions['electric']
        forcing=sum(g['electric'][i].T@(self.particles*physical_forcing[:,i]) for i in range(3))+thermal['reduced_forcing']
        mwdot=sum(g['electric_rate'][i].T@(self.particles*drift[:,i])+g['electric'][i].T@(self.particles*(g['electric_rate'][i]@w)) for i in range(3))
        rhs=forcing-self.mu/E**2*(g['ddot'].T@b)-self.mu/E*mwdot
        operator=h+self.mu/E**2*(d@solve(d.T.toarray()))
        bdot=positive_solve(operator,np.asarray(d@solve(rhs)));electric=solve(rhs-self.mu/E**2*(d.T@bdot))
        for iteration in range(4):
            r1=d@electric-h@bdot;r2=rhs-self.mu/E**2*(d.T@bdot)-mass@electric
            scale1=max(np.linalg.norm(d@electric)+np.linalg.norm(h@bdot),np.finfo(float).tiny)
            scale2=max(np.linalg.norm(rhs)+np.linalg.norm(self.mu/E**2*(d.T@bdot))+np.linalg.norm(mass@electric),np.finfo(float).tiny)
            if np.linalg.norm(r1)<2e-13*scale1 and np.linalg.norm(r2)<2e-13*scale2:break
            correction=positive_solve(operator,r1+d@solve(r2));bdot+=correction;electric+=solve(r2-self.mu/E**2*(d.T@correction))
        wdot=solve(-(d.T@bdot+g['ddot'].T@b)/E-mwdot)
        norm=lambda residual,*terms:float(np.linalg.norm(residual)/max(sum(np.linalg.norm(v) for v in terms),np.finfo(float).tiny))
        errors=dict(faraday=norm(h@bdot-d@electric,h@bdot,d@electric),ohm=norm(mass@electric-forcing-self.mu/E*(mass@wdot),mass@electric,forcing,self.mu/E*(mass@wdot)))
        if max(errors.values())>1e-10 or not np.isfinite(bdot).all():raise RuntimeError(f'Full vector material equations failed: {errors}')

        fields=g['fields'];volume=fields['volume'];bv=magnetic.reshape(self.shape+(3,));wv=drift.reshape(self.shape+(3,));r=fields['points'][...,0]
        directional_b=np.einsum('cqa,cqna->cqn',bv[...,[0,2]],fields['gradient'])
        local=np.einsum('cq,cqna->cna',volume*np.sum(bv*bv,axis=-1)/(2*MU0),fields['divergence'])
        local-=np.einsum('cq,cqa,cqn->cna',volume/MU0,bv[...,[0,2]],directional_b)
        local[...,0]-=np.einsum('cq,qn->cn',volume*bv[...,1]**2/(MU0*r),f.shape)
        local_torque=-np.einsum('cq,cqn->cn',volume*r*bv[...,1]/MU0,directional_b)
        directional_w=np.einsum('cqa,cqna->cqn',wv[...,[0,2]],fields['gradient'])
        local+=self.mu*np.einsum('cq,cqa,cqn->cna',f.particles,wv[...,[0,2]],directional_w)
        local[...,0]+=self.mu*np.einsum('cq,qn->cn',f.particles*wv[...,1]**2/r,f.shape)
        local_torque+=self.mu*np.einsum('cq,cqn->cn',f.particles*r*wv[...,1],directional_w)
        force=np.zeros_like(state.fluid.coordinates);np.add.at(force,f.ids.ravel(),local.reshape(-1,2));force[f.domain.axis,0]=0
        torque=np.zeros(f.nodes);np.add.at(torque,f.ids.ravel(),local_torque.ravel())

        joule=collisions['drag_power'];heating=collisions['heating']
        exchange=collisions['thermal_exchange_to_electrons'];frequency=collisions['exchange_frequency'];log=collisions['coulomb_log']
        entropy_rate=thermal['entropy_rate']+self.species.heating_rate(g,heating)
        internal=fields['specific_internal'].reshape(-1,2)
        thermal_power=float(np.sum(self.particles[:,None]*internal*entropy_rate.reshape(-1,2)))
        magnetic_power=float(b@h@bdot+.5*b@g['hdot']@b)
        current_power=float(self.mu*w@(mass@wdot)+.5*self.mu*w@mwdot)
        mechanical_power=float(np.sum(force*fields['velocity'])+torque@fields['omega'])
        kinetic_potential=self.species.potential_gradient(-self.delta*self.mu/(2*E)*np.sum(drift*drift,axis=-1))
        material_electric=electric+thermal['potential_electric']+kinetic_potential
        physical_electric=np.column_stack([v@material_electric for v in g['electric']])-np.cross(g['velocity'],magnetic)
        charge=self.space.gradient.T@(mass@w)
        result.update(magnetic_rate=bdot,electric_reduced=electric,electric_material=material_electric,electric_field=physical_electric,
            drift_rate=wdot,entropy_rate=entropy_rate,force=force,torque=torque,residuals=errors,
            magnetic_power=magnetic_power,current_power=current_power,mechanical_power=mechanical_power,thermal_relative_power=thermal_power,
            joule_power=float(g['volume']@joule),exchange_power_to_electrons=float(g['volume']@exchange),
            collision_electric=collisions['electric'],collision_heat_species_power=g['volume']@heating,
            collision_power_residual=float(g['volume']@(heating.sum(axis=1)-E*n*np.sum(drift*collisions['electric'],axis=1))),
            collision_entropy_rate_over_kB=float(np.sum(g['volume'][:,None]*heating/(E*fields['temperature'].reshape(-1,2)))),
            maximum_collision_relaxation_s=float(collisions['relaxation_rate'].max()),
            entropy_rates_over_kB=np.sum(self.particles[:,None]*entropy_rate.reshape(-1,2)/(f.gammas-1),axis=0),
            maximum_transport_outflow_s=thermal['maximum_transport_outflow_s'],
            exchange_frequency_range_s=None if frequency is None else [float(frequency.min()),float(frequency.max())],
            exchange_coulomb_log_range=None if log is None else [float(log.min()),float(log.max())],
            pair_power_residual=float(magnetic_power+current_power+mechanical_power+thermal_power),
            weak_charge_residual=float(np.linalg.norm(charge)/max(np.linalg.norm(self.space.gradient)*np.linalg.norm(mass@w),np.finfo(float).tiny)))
        if not np.isfinite(entropy_rate).all():raise ValueError('Nonfinite physical species evolution; no state repair')
        return result

    def derivative(self,state):
        response=self.response(state);fluid=self.fluid.derivative(state.fluid)
        fluid.momentum+=response['force'];fluid.angular_momentum+=response['torque'];fluid.entropy+=response['entropy_rate']
        return MaterialCurrentState(fluid,response['magnetic_rate'])
