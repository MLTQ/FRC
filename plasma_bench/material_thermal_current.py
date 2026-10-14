"""Moving finite-mass current with species entropy, resistance and heat exchange."""
import numpy as np
from scipy.constants import elementary_charge as E
from .material_carrier_current import MaterialCarrierCurrent, MaterialCurrentState, positive_solve
from .material_current_geometry import meridional_geometry
from .material_entropy_pair import MaterialEntropyPair
from .material_dual_entropy import MaterialDualEntropy
from .coulomb_exchange_wall_plasma import equilibration_rate


class MaterialThermalCurrent(MaterialCarrierCurrent):
    """General species entropies in the existing meridional/toroidal sector.

    eta is an explicitly prescribed scalar resistivity, not an inferred closure.
    Optional Maxwellian exchange is conditional, with no logarithm floor.
    """
    def __init__(self,fluid,space,*,resistivity=0.,coulomb_exchange=True,electron_heat_fraction=None,transport='upwind'):
        super().__init__(fluid,space)
        if not np.isfinite(resistivity) or resistivity<0:raise ValueError('Finite nonnegative prescribed resistivity required')
        self.resistivity=float(resistivity);self.coulomb_exchange=bool(coulomb_exchange)
        self.electron_heat_fraction=fluid.ion_mass/fluid.mass if electron_heat_fraction is None else float(electron_heat_fraction)
        if not 0<=self.electron_heat_fraction<=1:raise ValueError('Electron friction-heat fraction must be in [0,1]')
        if transport not in ('skew','upwind'):raise ValueError('Explicit skew or dual-volume upwind entropy transport required')
        self.transport=transport;self.species=(MaterialDualEntropy if transport=='upwind' else MaterialEntropyPair)(fluid,space)

    def geometry(self,state,fields=None):
        return meridional_geometry(self,state,fields)

    def response(self,state,*,rate=True):
        result=super().response(state,rate=rate)
        if not rate:return result
        thermal=self.species.pair(state.fluid,result,result['drift'])
        n=result['fields']['n'].ravel();wr,wz=result['drift'].T
        electric_resistance=self.resistivity*E*n[:,None]*result['drift']
        resistive_forcing=result['er'].T@(self.particles*electric_resistance[:,0])+result['ez'].T@(self.particles*electric_resistance[:,1])
        joule=self.resistivity*(E*n)**2*(wr*wr+wz*wz)
        temperature=result['fields']['temperature'].reshape(-1,2);te,ti=temperature.T
        exchange=np.zeros_like(n);frequency=np.zeros_like(n);log=None
        if self.coulomb_exchange:
            frequency,log=equilibration_rate(n,te,ti,self.fluid.ion_mass)
            exchange=1.5*E*n*frequency*(ti-te)
        heating=np.column_stack((self.electron_heat_fraction*joule+exchange,(1-self.electron_heat_fraction)*joule-exchange))
        internal=result['fields']['specific_internal'].reshape(-1,2)
        entropy_rate=thermal['entropy_rate']+self.species.heating_rate(result,heating)

        # Complete original mixed rows, now including pressure and resistance.
        # Starting from the parent solve is an algebraic initial guess only.
        h=result['h'];d=result['d'];mass=result['mass'];solve=result['solve'];w=result['drift_coefficients'];b=state.magnetic
        coefficient=self.delta*(result['magnetic_field']+self.mu/E*(result['curl']@w))-self.mu/E*(result['gradient'][:,0,1]-result['gradient'][:,1,0])
        forcing=result['er'].T@(-self.particles*wz*coefficient)+result['ez'].T@(self.particles*wr*coefficient)
        forcing+=thermal['reduced_forcing']+resistive_forcing
        mwdot=result['grdot'].T@(self.particles*wr)+result['gzdot'].T@(self.particles*wz)
        mwdot+=result['er'].T@(self.particles*(result['grdot']@w))+result['ez'].T@(self.particles*(result['gzdot']@w))
        rhs=forcing-self.mu/E**2*(result['ddot'].T@b)-self.mu/E*mwdot
        operator=h+self.mu/E**2*(d@solve(d.T.toarray()))
        bdot=result['magnetic_rate'].copy();electric=result['electric_reduced'].copy()
        for iteration in range(5):
            r1=d@electric-h@bdot;r2=rhs-self.mu/E**2*(d.T@bdot)-mass@electric
            scale1=max(np.linalg.norm(d@electric)+np.linalg.norm(h@bdot),np.finfo(float).tiny)
            scale2=max(np.linalg.norm(rhs)+np.linalg.norm(self.mu/E**2*(d.T@bdot))+np.linalg.norm(mass@electric),np.finfo(float).tiny)
            if np.linalg.norm(r1)<2e-13*scale1 and np.linalg.norm(r2)<2e-13*scale2:break
            correction=positive_solve(operator,r1+d@solve(r2));bdot+=correction
            electric+=solve(r2-self.mu/E**2*(d.T@correction))
        wdot=solve(-(d.T@bdot+result['ddot'].T@b)/E-mwdot)
        norm=lambda residual,*terms:float(np.linalg.norm(residual)/max(sum(np.linalg.norm(v) for v in terms),np.finfo(float).tiny))
        errors=dict(faraday=norm(h@bdot-d@electric,h@bdot,d@electric),
            ohm=norm(mass@electric-forcing-self.mu/E*(mass@wdot),mass@electric,forcing,self.mu/E*(mass@wdot)))
        if max(errors.values())>1e-10 or not np.isfinite(entropy_rate).all():
            raise RuntimeError(f'Physical thermal material equations failed: {errors}')
        magnetic_power=float(b@h@bdot+.5*b@result['hdot']@b)
        current_power=float(self.mu*w@(mass@wdot)+.5*self.mu*w@mwdot)
        thermal_power=float(np.sum(self.particles[:,None]*internal*entropy_rate.reshape(-1,2)))
        kinetic_potential=self.species.potential_gradient(-self.delta*self.mu/(2*E)*(wr*wr+wz*wz))
        material_electric=electric+thermal['potential_electric']+kinetic_potential
        velocity=result['fields']['sampled_velocity'].reshape(-1,2);magnetic=result['magnetic_field']
        physical_electric=np.column_stack((result['er']@material_electric+velocity[:,1]*magnetic,
            result['ez']@material_electric-velocity[:,0]*magnetic))
        result.update(magnetic_rate=bdot,electric_reduced=electric,electric_material=material_electric,electric_field=physical_electric,
            drift_rate=wdot,entropy_rate=entropy_rate,residuals=errors,magnetic_power=magnetic_power,current_power=current_power,
            thermal_relative_power=thermal_power,joule_power=float(result['volume']@joule),
            exchange_power_to_electrons=float(result['volume']@exchange),
            entropy_rates_over_kB=np.sum(self.particles[:,None]*entropy_rate.reshape(-1,2)/(self.fluid.gammas-1),axis=0),
            exchange_frequency_range_s=[float(frequency.min()),float(frequency.max())],
            exchange_coulomb_log_range=None if log is None else [float(log.min()),float(log.max())],
            maximum_transport_outflow_s=thermal.get('maximum_transport_outflow_s'),
            pair_power_residual=float(magnetic_power+current_power+result['mechanical_power']+thermal_power))
        return result

    def derivative(self,state):
        response=self.response(state);fluid=self.fluid.derivative(state.fluid)
        fluid.momentum+=response['force'];fluid.entropy+=response['entropy_rate']
        return MaterialCurrentState(fluid,response['magnetic_rate'])
