"""Prescribed finite coils driving the freely moving finite-mass plasma."""
import numpy as np
from scipy.constants import elementary_charge as E
from .electromagnetics import loop_field
from .material_vacuum_current import MaterialVacuumCurrent
from .material_carrier_current import MaterialCurrentState,positive_solve


class MaterialCoilCurrent(MaterialVacuumCurrent):
    def __init__(self,fluid,space,vacuum,coils,waveform,**options):
        super().__init__(fluid,space,vacuum,**options)
        self.coils=tuple(coils);self.waveform=waveform

    def source(self,coordinates,points):
        ids=self.fluid.domain.ids
        polygon=coordinates[np.r_[ids[:,0],ids[-1,1:],ids[-2::-1,-1]]]
        psi=np.zeros(len(points));field=np.zeros((len(points),3))
        for radius,center,turns in self.coils:
            # A prescribed filament may not cross into the evolving material.
            inside=False
            for a,b in zip(polygon,np.roll(polygon,-1,axis=0)):
                edge=b-a;parameter=np.clip(np.dot(np.array([radius,center])-a,edge)/np.dot(edge,edge),0.,1.)
                if np.linalg.norm(np.array([radius,center])-a-parameter*edge)<1e-12:
                    raise ValueError('Prescribed source filament touches the moving plasma')
                if (a[1]>center)!=(b[1]>center):
                    crossing=a[0]+(center-a[1])*(b[0]-a[0])/(b[1]-a[1])
                    if radius<crossing:inside=not inside
            if inside:raise ValueError('Prescribed source filament entered the moving plasma')
            p,br,bz=loop_field(points[:,0],points[:,1],radius,center)
            psi+=turns*p;field[:,0]+=turns*br;field[:,2]+=turns*bz
        potential=np.zeros_like(field);potential[:,1]=psi/points[:,0]
        return potential,field

    def response(self,state,*,rate=True,time=0.):
        # Retain the parent's actual thermal, force and geometry calculation.
        # Re-solve only the mixed electromagnetic rows with the exact source.
        a=super().response(state,rate=rate);f=self.fluid
        amplitude=float(self.waveform.value(time));slope=float(self.waveform.derivative(time))
        if not np.isfinite([amplitude,slope]).all():raise ValueError('Finite prescribed source waveform required')
        potential,unit=self.source(state.fluid.coordinates,a['fields']['points'].reshape(-1,2))
        maps=a['electric'];load=sum(g.T@(self.particles*potential[:,i]) for i,g in enumerate(maps))
        w=a['drift_coefficients'];drift=a['drift'];induced=a['magnetic_field'].copy();total=induced+amplitude*unit
        linkage=float(E*w@load)
        a.update(induced_magnetic_field=induced,magnetic_field=total,source_current_a=amplitude,source_current_rate_a_s=slope,
            source_linkage_wb=linkage,mutual_energy=amplitude*linkage,response_energy=a['total_energy']+amplitude*linkage,
            external_power=-slope*linkage)
        if not rate:return a
        h,d,mass,solve=a['h'],a['d'],a['mass'],a['solve'];b=state.magnetic
        thermal=self.species.pair(state.fluid,a,drift)
        curlw=np.column_stack([v@w for v in a['curl']])
        forcing_physical=self.delta*np.cross(drift,total+self.mu/E*curlw)-self.mu/E*np.cross(drift,a['vorticity'])
        n=a['fields']['n'].ravel();forcing_physical+=self.resistivity*E*n[:,None]*drift
        forcing=sum(maps[i].T@(self.particles*forcing_physical[:,i]) for i in range(3))+thermal['reduced_forcing']
        source_material=-slope*potential+amplitude*np.cross(a['velocity'],unit)
        drive=sum(maps[i].T@(self.particles*source_material[:,i]) for i in range(3))
        mwdot=sum(a['electric_rate'][i].T@(self.particles*drift[:,i])+maps[i].T@(self.particles*(a['electric_rate'][i]@w)) for i in range(3))
        rhs=forcing-drive-self.mu/E**2*(a['ddot'].T@b)-self.mu/E*mwdot
        operator=h+self.mu/E**2*(d@solve(d.T.toarray()))
        bdot=positive_solve(operator,np.asarray(d@solve(rhs)));electric=solve(rhs-self.mu/E**2*(d.T@bdot))
        for _ in range(4):
            r1=d@electric-h@bdot;r2=rhs-self.mu/E**2*(d.T@bdot)-mass@electric
            scale1=max(np.linalg.norm(d@electric)+np.linalg.norm(h@bdot),np.finfo(float).tiny)
            scale2=max(np.linalg.norm(rhs)+np.linalg.norm(self.mu/E**2*(d.T@bdot))+np.linalg.norm(mass@electric),np.finfo(float).tiny)
            if np.linalg.norm(r1)<2e-13*scale1 and np.linalg.norm(r2)<2e-13*scale2:break
            correction=positive_solve(operator,r1+d@solve(r2));bdot+=correction;electric+=solve(r2-self.mu/E**2*(d.T@correction))
        wdot=solve(-(d.T@bdot+a['ddot'].T@b)/E-mwdot);projected_material=electric+solve(drive)
        norm=lambda r,*v:float(np.linalg.norm(r)/max(sum(np.linalg.norm(t) for t in v),np.finfo(float).tiny))
        errors=dict(faraday=norm(h@bdot-d@electric,h@bdot,d@electric),ohm=norm(mass@projected_material-forcing-self.mu/E*(mass@wdot),mass@projected_material,forcing,self.mu/E*(mass@wdot)))
        if max(errors.values())>1e-10 or not np.isfinite(bdot).all():raise RuntimeError(f'Driven moving plasma equations failed: {errors}')
        force_density=E*n[:,None]*np.cross(drift,amplitude*unit);vol=a['fields']['volume']
        force_density=force_density.reshape(self.shape+(3,));r=a['fields']['points'][...,0]
        local=np.einsum('cq,cqa,qn->cna',vol,force_density[...,[0,2]],f.shape)
        torque=np.einsum('cq,qn->cn',vol*r*force_density[...,1],f.shape)
        force=np.zeros_like(state.fluid.coordinates);spin=np.zeros(f.nodes)
        np.add.at(force,f.ids.ravel(),local.reshape(-1,2));np.add.at(spin,f.ids.ravel(),torque.ravel());force[f.domain.axis,0]=0
        mechanical=a['mechanical_power']+float(np.sum(force*a['fields']['velocity'])+spin@a['fields']['omega'])
        magnetic=float(b@h@bdot+.5*b@a['hdot']@b);current=float(self.mu*w@(mass@wdot)+.5*self.mu*w@mwdot)
        potential_rate=np.zeros_like(potential)
        potential_rate[:,1]=a['velocity'][:,0]*(unit[:,2]-potential[:,1]/a['r'])-a['velocity'][:,2]*unit[:,0]
        load_rate=sum(a['electric_rate'][i].T@(self.particles*potential[:,i])+maps[i].T@(self.particles*potential_rate[:,i]) for i in range(3))
        linkage_rate=float(E*(wdot@load+w@load_rate))
        delta_e=electric-a['electric_reduced']
        a['electric_field']+=np.column_stack([g@delta_e for g in maps])-slope*potential
        a['force']+=force;a['torque']+=spin
        a.update(magnetic_rate=bdot,drift_rate=wdot,electric_reduced=electric,electric_induced_material=a['electric_material']+delta_e,
            electric_material=a['electric_material']+delta_e+solve(drive),residuals=errors,magnetic_power=magnetic,current_power=current,mechanical_power=mechanical,
            source_linkage_rate_wb_s=linkage_rate,source_response_power=amplitude*linkage_rate,
            pair_power_residual=magnetic+current+mechanical+a['thermal_relative_power']+slope*linkage)
        return a

    def derivative(self,state,time=0.):
        a=self.response(state,time=time);fluid=self.fluid.derivative(state.fluid)
        fluid.momentum+=a['force'];fluid.angular_momentum+=a['torque'];fluid.entropy+=a['entropy_rate']
        return MaterialCurrentState(fluid,a['magnetic_rate'])

    def step(self,state,dt,time=0.):
        if not np.isfinite(dt) or dt<=0:raise ValueError('Positive finite step required')
        k1=self.derivative(state,time);k2=self.derivative(state.added(k1,.5*dt),time+.5*dt)
        k3=self.derivative(state.added(k2,.5*dt),time+.5*dt);k4=self.derivative(state.added(k3,dt),time+dt)
        result=state.added(k1,dt/6).added(k2,dt/3).added(k3,dt/3).added(k4,dt/6)
        self.response(result,rate=False,time=time+dt);return result
