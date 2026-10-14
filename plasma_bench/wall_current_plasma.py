"""Moving collisional plasma with a physical evolved toroidal wall current."""
from copy import copy
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME, mu_0 as MU0
from .regular_plasma import RegularPlasma, RegularPlasmaRate
from .regular_species import RegularSpecies
from .split_inverse_wall_induction import SplitInverseWallInduction
from .wall_current_geometry import WallCurrentTransport
from .wall_current_model import WallCurrentModel
from .variable_compatible_current import weak_vector
from .coulomb_exchange_wall_plasma import equilibration_rate
from .spitzer_wall_plasma import spitzer_samples


@dataclass(frozen=True)
class WallCurrentRate(RegularPlasmaRate):
    wall_current: np.ndarray


class WallCurrentPlasma(RegularPlasma):
    def __init__(self, mapping, ion_mass, *, drive, transport='fixed', resistivity=None, **options):
        super().__init__(mapping,ion_mass,**options)
        if transport not in ('fixed','spitzer'):
            raise ValueError('Explicit fixed or spitzer transport required')
        if transport=='fixed' and (resistivity is None or not np.isfinite(resistivity) or resistivity<0):
            raise ValueError('Nonnegative prescribed fixed resistance required')
        if not np.allclose(self.fluid.thermal.gammas,5/3,rtol=0.,atol=4*np.finfo(float).eps):
            raise ValueError('Both physical species require gamma=5/3')
        self.drive=drive;self.transport=transport;self.resistivity=resistivity

    def field_value(self, when):return float(self.drive.value(when))

    def operators(self, model, reduced):
        if self._operators is None or not self.cache_geometry:
            items=(WallCurrentTransport(model),RegularSpecies(model),
                SplitInverseWallInduction(reduced,**self.solver_options))
            if not self.cache_geometry:return items
            self._operators=items
        a,b,c=(copy(x) for x in self._operators);a.model=b.model=model;c.model=reduced
        return a,b,c

    def derivative(self, state, field, when):
        state=np.asarray(state,float);g=self.geometry
        if not np.isfinite(when) or field.axial_bias!=self.field_value(when):
            raise ValueError('Actual field must match the explicit source time')
        fields=self.fluid.fields(state);packed=field.packed_increment(g)
        model=WallCurrentModel(self.mapping,state[...,0],self.ion_mass,field.wall_current)
        # Maps/preconditioner are cached; the actual Schur inertia is rebuilt below.
        transport_op,species_op,induction=self.operators(model,model.base)
        bulk,_=self.fluid.derivative(state);nt=bulk[...,0]
        electromagnetic=model.electromagnetic_pair(state,field,nt);current=electromagnetic['current']
        transport=transport_op.pair(state,current,nt);species=species_op.pair(state,current)
        j=g.sample_electric(current);b=field.samples(g)
        if self.transport=='spitzer':
            n,temperature,eta,frequency,log=spitzer_samples(self.mapping,state,self.ion_mass)
        else:
            eta=np.full_like(model.density,self.resistivity)
        heat=eta*np.sum(j*j,axis=-1)
        fluid=bulk+electromagnetic['fluid_rate']+transport['fluid_rate']+species['fluid_rate']
        fluid[...,4]+=self.mapping.adjoint(heat)/(self.fluid.volume*fields['specific_internal'][...,0])
        te,ti=np.moveaxis(fields['temperature'],-1,0)
        exchange_frequency,exchange_log=equilibration_rate(fields['n'],te,ti,self.ion_mass)
        q=1.5*E*fields['n']*exchange_frequency*(ti-te)
        fluid[...,4:]+=np.stack((q,-q),axis=-1)/fields['specific_internal']
        reaction=-fluid[-1,-1,:,1].copy();fluid[-1,-1,:,1]=0.
        source_rate=float(self.drive.derivative(when))
        external=np.zeros_like(j);external[...,1]=-.5*g.r*source_rate
        extra=weak_vector(g,eta[...,None]*j-external)
        forcing=electromagnetic['forcing']+transport['forcing']+species['reduced_forcing']+extra
        reduced,right,reconstruct=model.reduced_induction(forcing)
        induction.model=reduced;solved=induction.solve(right)
        current_rate,wall_rate=reconstruct(solved.current)
        electric=g.embed(solved.electric)+species['potential_electric']
        full_forcing=electromagnetic['forcing']+transport['forcing']+species['forcing']+extra
        inertia=model.inertia@current_rate;emass=model.mass@electric
        scale=np.linalg.norm(abs(emass)+abs(inertia)+abs(full_forcing))
        backward=float(np.linalg.norm(emass-inertia-full_forcing)/max(scale,np.finfo(float).tiny))
        residual=max(solved.residual,backward)
        if not np.isfinite(residual) or residual>self.solver_options['tolerance']:
            raise RuntimeError(f'Full wall-current Ohm residual {backward} failed; no input modified')
        if not all(np.isfinite(x).all() for x in (fluid,electric,solved.magnetic,current_rate,wall_rate)):
            raise ValueError('Nonfinite physical wall-current rate')
        grid_fluid=float(np.sum(self.fluid.volume[...,None]*fields['energy_gradient']*fluid))
        current_power=float(np.sum(model.mesh.volumes*model.current_energy_rate(current,current_rate,nt)['total']))
        harmonic=float(np.sum(model.mesh.volumes)*field.axial_bias*source_rate/MU0)
        magnetic_power=float(np.sum(model.mesh.volumes*model.magnetic_energy_rate(field,solved.magnetic)))+harmonic
        powers=dict(fluid=grid_fluid,current=current_power,magnetic=magnetic_power)
        source_plasma=float(np.sum(g.weight*np.sum(j*external,axis=-1)));source=source_plasma+harmonic
        diagnostics=dict(powers=powers,total_power=sum(powers.values()),source_power_w=source,
            source_plasma_power_w=source_plasma,source_harmonic_power_w=harmonic,
            power_balance_residual_w=sum(powers.values())-source,full_species_ohm_backward=backward,
            induction_residuals=solved.residual_components,iterations=solved.iterations,
            wall_normal_reaction_max=float(abs(reaction).max()),
            wall_current_range_a_m2=[float(field.wall_current.min()),float(field.wall_current.max())],
            wall_loop_voltage_v=float(-np.pi*model.mesh.radius**2*source_rate),
            ohmic_power=float(np.sum(g.weight*heat)),resistivity_range_ohm_m=[float(eta.min()),float(eta.max())],
            exchange_power_to_electrons=float(np.sum(self.fluid.volume*q)),
            exchange_frequency_range_s=[float(exchange_frequency.min()),float(exchange_frequency.max())],
            exchange_coulomb_log_range=[float(exchange_log.min()),float(exchange_log.max())],
            electric_representation='E_ind=E_physical-E_ext; add exact E_ext,theta=-r*Bbar_dot/2. E_ind wall trace is zero; physical current wall trace evolves.')
        if self.transport=='fixed':diagnostics['resistivity_ohm_m']=self.resistivity
        else:
            hall=E*np.linalg.norm(b,axis=-1)/(ME*frequency)
            drift=np.linalg.norm(j,axis=-1)/(E*n)
            ratio=drift[...,None]/np.sqrt(E*temperature/np.array([ME,self.ion_mass]))
            diagnostics.update(maximum_electron_hall_parameter=float(hall.max()),
                maximum_drift_over_species_thermal_speed=ratio.reshape(-1,2).max(axis=0).tolist(),
                weak_field_review=bool(np.any(hall>=.1)),relative_drift_review=bool(np.any(ratio>=.1)))
        return WallCurrentRate(fluid,solved.magnetic,current_rate,electric,residual,diagnostics,wall_rate)
