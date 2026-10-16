"""Complete warm neutral compression and rebound against an energy integral."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, elementary_charge as EV
from scipy.integrate import quad
from plasma_bench.species_material import MaterialSpecies
from plasma_bench.mhd_benchmarks import sod_exact
from .implicit_plasma_records import provenance, freeze


def run(destination):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh compression output required')
    identity=provenance(['examples.species_material_compression']);out.mkdir(parents=True);freeze(out,identity)
    start=time.monotonic();a=.025;density=1e19;temperature=1.;mass=39.948*atomic_mass;gamma=5/3;smin=.5
    acceleration=2*EV*temperature/(mass*a*a)
    def speed(scale):
        return np.sqrt(2*acceleration/(gamma-1)*smin**(1-gamma)*(-np.expm1((1-gamma)*np.log(scale/smin))))
    def integrand(y):
        return 2*y/speed(smin+y*y) if y else np.sqrt(2/(acceleration*smin**-gamma))
    turn=quad(integrand,0.,np.sqrt(1-smin),epsabs=1e-15,epsrel=1e-11)[0]
    rate=speed(1.);duration=2*turn
    report=dict(terminal=False,source_sha256=identity['source_sha256'],duration_s=duration,turn_time_s=turn,
                minimum_scale=smin,initial_scale_rate_s=-rate,cases=[],
                scope='Neutral argon Euler compression/rebound with quadratic pressure and exact vacuum. Independent scalar energy-integral duration, no forced trajectory; numerical viscosity comparison, not lab validation.')
    def save():
        report['elapsed_seconds']=time.monotonic()-start;(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save();fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    refscale=np.linspace(smin,1.,101)
    elapsed=np.array([quad(integrand,0.,np.sqrt(s-smin),epsabs=1e-15,epsrel=1e-10)[0] for s in refscale])
    axes[0].plot(np.r_[turn-elapsed[::-1],turn+elapsed[1:]]*1e6,np.r_[refscale[::-1],refscale[1:]],'k--',label='Continuum energy integral')
    for limiter in ('none','minmod'):
        for cells in (32,64,128):
            faces=-a*np.cos(np.linspace(0,np.pi,cells+1))[None,:];numbers=density*np.diff(faces)
            temp=temperature*(1-(faces[:,:-1]**2+faces[:,:-1]*faces[:,1:]+faces[:,1:]**2)/(3*a*a))
            velocity=np.zeros((1,cells+1,3));velocity[...,2]=-rate*faces
            model=MaterialSpecies(numbers,[mass],[0.],(-.2,.2),numerical_viscosity=(.5,2),viscosity_limiter=limiter)
            state=model.state(faces,temp,velocity);initial=model.totals(state);times=[0.];states=[state];energies=[initial['energy_j_m2']];heat=0.;failure=None
            def power(s):return float(np.sum(-model.pressure(s)[1]*np.diff(s.velocity[...,2],axis=1)))
            while times[-1]<duration:
                dt=min(model.timestep(state),duration-times[-1])
                try:candidate=model.step(state,dt)
                except (ValueError,RuntimeError) as error:
                    failure=dict(time_s=times[-1],dt_s=dt,message=str(error));break
                heat+=dt/2*(power(state)+power(candidate));state=candidate
                times.append(min(times[-1]+dt,duration));states.append(state);energies.append(model.totals(state)['energy_j_m2'])
            name=f'{limiter}_n{cells}';energy=np.asarray(energies)
            record=dict(name=name,limiter=limiter,cells=cells,success=failure is None,failure=failure,steps=len(times)-1,completed_time_s=times[-1],
                maximum_position_error_m=float(abs(state.faces-faces).max()),
                velocity_relative_l2=float(np.linalg.norm(state.velocity+velocity)/np.linalg.norm(velocity)),
                thermal_relative_l1=float(abs(state.thermal-states[0].thermal).sum()/states[0].thermal.sum()),
                numerical_heat_j_m2=heat,maximum_relative_energy_error=float(abs(energy-energy[0]).max()/energy[0]),
                minimum_thermal_j_m2=float(min(s.thermal.min() for s in states)))
            report['cases'].append(record)
            np.savez_compressed(out/(name+'.npz'),time_s=times,faces_m=[s.faces for s in states],velocity_m_s=[s.velocity for s in states],
                                thermal_j_m2_history=[s.thermal for s in states],energy_j_m2=energy,numbers_m2=numbers)
            axes[0].plot(np.asarray(times)*1e6,[s.faces[0,-1]/a for s in states],label=f'{limiter}, {cells}',alpha=.75)
            axes[1].plot(np.asarray(times)*1e6,(energy-energy[0])/energy[0],label=f'{limiter}, {cells}')
            save();print(record,flush=True)
    axes[0].set(xlabel='Time [µs]',ylabel='Physical front / initial front',title='Full compression and rebound')
    axes[1].set(xlabel='Time [µs]',ylabel='Total energy change / initial',title='Actual energy conservation')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.savefig(out/'compression.png',dpi=170);fig.savefig(out/'compression.pdf');plt.close(fig)
    if identity!=provenance(['examples.species_material_compression']):raise RuntimeError('Source changed during compression run')
    report['terminal']=True;report['success']=all(c['success'] for c in report['cases']);save();return report


def run_shock(destination):
    """Finite neutral slab; compare its central shock before edge fans arrive."""
    out=Path(destination)
    if out.exists():raise ValueError('Fresh shock output required')
    identity=provenance(['examples.species_material_compression']);out.mkdir(parents=True);freeze(out,identity)
    start=time.monotonic();a=.05;density=1e19;mass=39.948*atomic_mass;speed0=np.sqrt(EV/mass);tau=.12;duration=tau*a/speed0;gamma=5/3
    report=dict(terminal=False,source_sha256=identity['source_sha256'],duration_s=duration,cases=[],
                scope='Finite free neutral argon slab with a central Sod discontinuity. Compare only |x|<a/2 at t=.12a/sqrt(eV/m), before free-edge rarefactions reach this region; exact vacuum and full-domain energy retained.')
    def save():
        report['elapsed_seconds']=time.monotonic()-start;(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save();fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    xplot=np.linspace(-.5,.5,4001);exact=sod_exact(xplot,tau,gamma)
    for ax,values in zip(axes,exact):ax.plot(xplot,values,'k--',label='Exact Euler Riemann solution')
    _,vstar,pstar=sod_exact(np.array([0.]),tau,gamma);vstar=float(vstar[0]);pstar=float(pstar[0])
    shock=np.sqrt(gamma*.1/.125)*np.sqrt((gamma+1)/(2*gamma)*pstar/.1+(gamma-1)/(2*gamma))
    breaks=tau*np.array([-np.sqrt(gamma),vstar-np.sqrt(gamma)*pstar**((gamma-1)/(2*gamma)),vstar,shock])*a
    nodes,weights=np.polynomial.legendre.leggauss(16)
    for limiter in ('none','minmod'):
        for cells in (64,128,256):
            faces=np.linspace(-a,a,cells+1)[None,:];cent=(faces[:,:-1]+faces[:,1:])/2
            number=density*np.where(cent<0,1.,.125);temp=np.where(cent<0,1.,.8)
            model=MaterialSpecies(number*np.diff(faces),[mass],[0.],(-.4,.4),numerical_viscosity=(.5,2),viscosity_limiter=limiter)
            state=model.state(faces,temp);times=[0.];states=[state];initial=model.totals(state);energies=[initial['energy_j_m2']];momenta=[initial['momentum_kg_per_m_s']];heat=0.;failure=None
            def power(s):return float(np.sum(-model.pressure(s)[1]*np.diff(s.velocity[...,2],axis=1)))
            while times[-1]<duration:
                dt=min(model.timestep(state),duration-times[-1])
                try:candidate=model.step(state,dt)
                except (ValueError,RuntimeError) as error:
                    failure=dict(time_s=times[-1],dt_s=dt,message=str(error));break
                heat+=dt/2*(power(state)+power(candidate));state=candidate
                times.append(min(times[-1]+dt,duration));states.append(state)
                totals=model.totals(state);energies.append(totals['energy_j_m2']);momenta.append(totals['momentum_kg_per_m_s'])
            edges=np.unique(np.r_[-a/2,a/2,breaks,state.faces[0,(state.faces[0]>-a/2)&(state.faces[0]<a/2)]])
            x=((edges[:-1,None]+edges[1:,None])/2+np.diff(edges)[:,None]*nodes/2).ravel()
            w=(np.diff(edges)[:,None]*weights/2).ravel()/a;sample=model.sample(state,x)
            actual=[sample['number'][0]/density,sample['velocity'][0,:,2]/speed0,sample['number'][0]*sample['temperature_ev'][0]/density]
            reference=sod_exact(x/a,times[-1]*speed0/a,gamma)
            errors=[float(np.sum(w*abs(v-r))) for v,r in zip(actual,reference)]
            name=f'{limiter}_n{cells}';energy=np.asarray(energies)
            record=dict(name=name,limiter=limiter,cells=cells,success=failure is None,failure=failure,steps=len(times)-1,completed_time_s=times[-1],
                central_normalized_l1=dict(zip(('density','velocity','pressure'),errors)),numerical_heat_j_m2=heat,
                maximum_relative_energy_error=float(abs(energy-energy[0]).max()/energy[0]),
                maximum_momentum_change_kg_per_m_s=float(np.linalg.norm(np.asarray(momenta)-momenta[0],axis=1).max()),
                minimum_thermal_j_m2=float(min(s.thermal.min() for s in states)))
            np.savez_compressed(out/(name+'.npz'),time_s=times,faces_m=[s.faces for s in states],velocity_m_s=[s.velocity for s in states],
                                thermal_j_m2_history=[s.thermal for s in states],energy_j_m2=energy,momentum_kg_per_m_s=momenta,numbers_m2=model.numbers)
            if cells==256:
                sample=model.sample(state,xplot*a)
                for ax,values in zip(axes,[sample['number'][0]/density,sample['velocity'][0,:,2]/speed0,sample['number'][0]*sample['temperature_ev'][0]/density]):
                    ax.plot(xplot,values,label=f'{limiter}, {cells}')
            report['cases'].append(record);save();print(record,flush=True)
    for ax,title in zip(axes,('Density / left density','Velocity / sqrt(eV/m)','Pressure / left pressure')):
        ax.set(xlabel='Position / initial half-width',ylabel=title);ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.savefig(out/'shock.png',dpi=170);fig.savefig(out/'shock.pdf');plt.close(fig)
    if identity!=provenance(['examples.species_material_compression']):raise RuntimeError('Source changed during shock run')
    report['terminal']=True;report['success']=all(c['success'] for c in report['cases']);save();return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');parser.add_argument('--shock',action='store_true');args=parser.parse_args()
    result=(run_shock if args.shock else run)(args.destination);raise SystemExit(0 if result['success'] else 1)
