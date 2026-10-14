"""Moving finite-tube plasma driven by finite coils through an empty exterior."""
from copy import copy
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME, mu_0 as MU0
from .regular_plasma import RegularPlasmaRate
from .regular_current import RegularCurrent
from .finite_cylinder_map import FiniteSpecies
from .wall_current_plasma import WallCurrentPlasma
from .finite_current_transport import FiniteCurrentTransport
from .finite_induction import FiniteInduction
from .variable_compatible_current import weak_vector
from .coulomb_exchange_wall_plasma import equilibration_rate
from .spitzer_wall_plasma import spitzer_samples


class FiniteCoilPlasma(WallCurrentPlasma):
    def __init__(self,mapping,ion_mass,*,source,**options):
        super().__init__(mapping,ion_mass,**options)
        self.source=source

    def operators(self, model):
        if self._operators is None or not self.cache_geometry:
            items=(FiniteCurrentTransport(model),FiniteSpecies(model),FiniteInduction(model,**self.solver_options))
            if not self.cache_geometry:return items
            self._operators=items
        result=tuple(copy(x) for x in self._operators)
        for x in result:x.model=model
        return result

    def derivative(self, state, field, when):
        state=np.asarray(state,float); g=self.geometry
        if not np.isfinite(when) or field.source_amplitude != self.field_value(when):
            raise ValueError('Actual field must match the source time')
        fields=self.fluid.fields(state); packed=field.packed_increment(g)
        model=RegularCurrent(self.mapping,state[...,0],self.ion_mass)
        transport_op,species_op,induction=self.operators(model)
        bulk,_=self.fluid.derivative(state); nt=bulk[...,0]
        electromagnetic=model.electromagnetic_pair(state,field,nt); current=electromagnetic['current']
        transport=transport_op.pair(state,current,nt); species=species_op.pair(state,current)
        j=g.sample_electric(current); b=field.samples(g)
        if self.transport=='spitzer':
            n,temperature,eta,frequency,log=spitzer_samples(self.mapping,state,self.ion_mass)
        else:eta=np.full_like(model.density,self.resistivity)
        heat=eta*np.sum(j*j,axis=-1)
        fluid=bulk+electromagnetic['fluid_rate']+transport['fluid_rate']+species['fluid_rate']
        fluid[...,4]+=self.mapping.adjoint(heat)/(self.fluid.volume*fields['specific_internal'][...,0])
        te,ti=np.moveaxis(fields['temperature'],-1,0)
        exchange_frequency,exchange_log=equilibration_rate(fields['n'],te,ti,self.ion_mass)
        q=1.5*E*fields['n']*exchange_frequency*(ti-te)
        fluid[...,4:]+=np.stack((q,-q),axis=-1)/fields['specific_internal']
        reaction=-fluid[-1,-1,:,1].copy(); fluid[-1,-1,:,1]=0.
        end_reaction=-fluid[...,[0,-1],3].copy();fluid[...,[0,-1],3]=0.
        source_rate=float(self.drive.derivative(when))
        external=-source_rate*self.source.a_unit
        extra=weak_vector(g,eta[...,None]*j-external)
        forcing=electromagnetic['forcing']+transport['forcing']+species['reduced_forcing']+extra
        solved=induction.solve(forcing)
        electric=solved.electric+species['potential_electric']
        full_forcing=electromagnetic['forcing']+transport['forcing']+species['forcing']+extra
        inertia=model.inertia@solved.current; emass=model.mass@electric
        scale=np.linalg.norm(abs(emass)+abs(inertia)+abs(full_forcing))
        backward=float(np.linalg.norm(emass-inertia-full_forcing)/max(scale,np.finfo(float).tiny))
        residual=max(solved.residual,backward)
        if not np.isfinite(residual) or residual>self.solver_options['tolerance']:
            raise RuntimeError(f'Full vacuum-coupled Ohm residual {backward} failed; no input modified')
        if not all(np.isfinite(x).all() for x in (fluid,electric,solved.magnetic,solved.current)):
            raise ValueError('Nonfinite physical vacuum-plasma rate')
        grid_fluid=float(np.sum(self.fluid.volume[...,None]*fields['energy_gradient']*fluid))
        current_power=float(np.sum(model.mesh.volumes*model.current_energy_rate(current,solved.current,nt)['total']))
        flux=self.source.linkage(current);flux_rate=self.source.linkage(solved.current)
        cross_power=source_rate*flux+field.source_amplitude*flux_rate
        induced=float(packed@(g.hessian@solved.magnetic))
        powers=dict(fluid=grid_fluid,current=current_power,magnetic=induced+cross_power)
        source=field.source_amplitude*flux_rate
        source_plasma=float(np.sum(g.weight*np.sum(j*external,axis=-1)))
        diagnostics=dict(powers=powers,total_power=sum(powers.values()),source_power_w=float(source),
            source_plasma_power_w=source_plasma,source_background_power_w=None,
            source_cross_power_w=float(cross_power),power_balance_residual_w=float(sum(powers.values())-source),
            full_species_ohm_backward=backward,induction_residuals=solved.residual_components,iterations=solved.iterations,
            wall_normal_reaction_max=float(abs(reaction).max()),
            source_linked_flux_wb=flux,source_linked_flux_rate_v=flux_rate,
            end_normal_reaction_max=float(abs(end_reaction).max()),
            ohmic_power=float(np.sum(g.weight*heat)),resistivity_range_ohm_m=[float(eta.min()),float(eta.max())],
            exchange_power_to_electrons=float(np.sum(self.fluid.volume*q)),
            exchange_frequency_range_s=[float(exchange_frequency.min()),float(exchange_frequency.max())],
            exchange_coulomb_log_range=[float(exchange_log.min()),float(exchange_log.max())],
            electric_representation='E_ind=E_physical-E_source; E_source=-Idot*Aunit from finite circular coils. Side/end flux evolves; vacuum contains no matter. Energy/work subtract prescribed filament self baseline.')
        if self.transport=='fixed':diagnostics['resistivity_ohm_m']=self.resistivity
        else:
            hall=E*np.linalg.norm(b,axis=-1)/(ME*frequency); drift=np.linalg.norm(j,axis=-1)/(E*n)
            ratio=drift[...,None]/np.sqrt(E*temperature/np.array([ME,self.ion_mass]))
            diagnostics.update(maximum_electron_hall_parameter=float(hall.max()),
                maximum_drift_over_species_thermal_speed=ratio.reshape(-1,2).max(axis=0).tolist(),
                weak_field_review=bool(np.any(hall>=.1)),relative_drift_review=bool(np.any(ratio>=.1)))
        return RegularPlasmaRate(fluid,solved.magnetic,solved.current,electric,residual,diagnostics)
