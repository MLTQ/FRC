"""Neutral, charged and electron-ion free expansions on independent supports."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, electron_mass, elementary_charge as EV, epsilon_0 as EPS0
from plasma_bench.species_material import MaterialSpecies, SpeciesMaterialState
from .atomic_moving import charged_expansion_time
from .implicit_plasma_records import provenance, freeze


def reference_density(points, kind, time_s, half_width, sound, scale):
    if kind=='charged':return (abs(points)<half_width*scale).astype(float)/scale
    if time_s==0:return (abs(points)<half_width).astype(float)
    xi=(abs(points)-half_width)/time_s
    ratio=np.where(xi<=-sound,1.,np.maximum(0.,.75*(1-xi/(3*sound))))
    return ratio**3


def density_error(model,state,kind,time_s,a,sound,scale,density):
    front=a*scale if kind=='charged' else a+3*sound*time_s
    head=a if kind=='charged' else a-sound*time_s
    edges=np.unique(np.r_[state.faces[0],-front,front,-head,head])
    node,weight=np.polynomial.legendre.leggauss(16)
    x=(edges[:-1,None]+edges[1:,None])/2+np.diff(edges)[:,None]*node/2
    sampled=model.sample(state,x.ravel())['number'][0].reshape(x.shape)/density
    exact=reference_density(x,kind,time_s,a,sound,scale)
    return float(np.sum(np.diff(edges)[:,None]*weight/2*abs(sampled-exact))/(2*a))


def run(destination, *, mesh='uniform', viscosity=(0.,0.), step_scale=1., grids=(32,64,128), viscosity_limiter='none'):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh independent-support output required')
    if mesh not in ('uniform','cosine'):raise ValueError('Choose uniform or cosine material labels')
    if not np.isfinite(step_scale) or not 0<step_scale<=1:raise ValueError('Time-step scale must lie in (0,1]')
    if len(grids)<2 or list(grids)!=sorted(set(grids)) or min(grids)<2:raise ValueError('At least two increasing cell counts required')
    identity=provenance(['examples.species_material']);out.mkdir(parents=True);freeze(out,identity)
    began=time.monotonic(); report=dict(terminal=False,source_sha256=identity['source_sha256'],cases=[],mesh=mesh,step_scale=step_scale,numerical_viscosity=list(viscosity),viscosity_limiter=viscosity_limiter,grids=list(grids),
        scope='Planar independent-support Euler-Poisson fluids, physical masses, no reactions/collisions; optional declared compressive numerical viscosity. Actual fields and front convergence, not formation or kinetic validation.')
    def save():
        report['elapsed_seconds']=time.monotonic()-began;(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save();finals={}; neutral_mass=39.948*atomic_mass; ion_mass=neutral_mass-electron_mass
    for kind in ('neutral','charged','plasma'):
        a=.01 if kind=='plasma' else .025
        density={'neutral':1e19,'charged':1e12,'plasma':1e14}[kind]
        temperature=1. if kind=='charged' else .1
        masses=[ion_mass,electron_mass] if kind=='plasma' else [ion_mass if kind=='charged' else neutral_mass]
        charge=[1.,-1.] if kind=='plasma' else [1. if kind=='charged' else 0.]
        sound=np.sqrt((5/3)*EV*temperature/masses[0]);scale=1.5
        duration=(charged_expansion_time(scale,2*EV*temperature/(ion_mass*a*a),EV**2*density/(EPS0*ion_mass))[0]
                  if kind=='charged' else (20e-9 if kind=='plasma' else a/(4*sound)))
        for cells in grids:
            labels=np.linspace(-1.,1.,cells+1) if mesh=='uniform' else -np.cos(np.linspace(0,np.pi,cells+1))
            faces=np.broadcast_to(a*labels,(len(masses),cells+1)).copy()
            numbers=density*np.diff(faces,axis=1)
            model=MaterialSpecies(numbers,masses,charge,(-.2,.2),left_electric=-EV*density*a/EPS0 if kind=='charged' else 0.,
                                  cfl=.15,maximum_plasma_phase=.04*32/cells,numerical_viscosity=viscosity,viscosity_limiter=viscosity_limiter)
            if kind=='charged':
                left,right=faces[:,:-1],faces[:,1:];temp=temperature*(1-(left*left+left*right+right*right)/(3*a*a))
            elif kind=='plasma':temp=np.array([.2,2.])[:,None]
            else:temp=temperature
            state=model.state(faces,temp); t=0.;times=[0.];states=[state];failure=None
            initial=model.totals(state);rows=[initial];viscous_heat=0.
            def viscous_power(s):return float(np.sum(-model.pressure(s)[1]*np.diff(s.velocity[...,2],axis=1)))
            while t<duration:
                dt=min(step_scale*model.timestep(state),duration-t)
                try:candidate=model.step(state,dt)
                except (ValueError,RuntimeError) as error:
                    failure=dict(step=len(times),time_s=t,dt_s=dt,message=str(error));break
                viscous_heat+=dt/2*(viscous_power(state)+viscous_power(candidate))
                state=candidate;t=min(t+dt,duration);times.append(t);states.append(state);rows.append(model.totals(state))
            budgets={k:np.asarray([r[k] for r in rows]) for k in rows[0]}
            energy=budgets['energy_j_m2']; momentum=budgets['momentum_kg_per_m_s']
            case=dict(name=f'{kind}_n{cells}',kind=kind,cells=cells,species_charge=charge,species_mass_kg=masses,
                reference_density_m3=density,duration_s=duration,completed_time_s=t,steps=len(times)-1,success=failure is None,failure=failure,
                cfl=model.cfl,maximum_plasma_phase=model.maximum_plasma_phase,numerical_viscosity=model.viscosity.tolist(),viscosity_limiter=model.viscosity_limiter,
                final_fronts_m=state.faces[:,[0,-1]].tolist(),
                maximum_relative_energy_error=float(abs(energy-energy[0]).max()/energy[0]),
                maximum_momentum_change_kg_per_m_s=float(np.linalg.norm(momentum-momentum[0],axis=1).max()),
                numerical_viscous_heat_j_m2=viscous_heat,
                maximum_particle_number_change=float(abs(budgets['numbers']-budgets['numbers'][0]).max()),
                minimum_thermal_j_m2=float(min(s.thermal.min() for s in states)),
                final_peak_electric_v_m=float(abs(model.field(state)[1]).max()))
            if kind!='plasma' and failure is None:
                front=a*scale if kind=='charged' else a+3*sound*t
                case.update(analytic_front_m=front,maximum_front_error_m=float(abs(state.faces[0,[0,-1]]-[-front,front]).max()),
                    final_number_l1_over_initial_number=density_error(model,state,kind,t,a,sound,scale,density))
                if kind=='charged':
                    edges,e,_,_=model.field(state);exact=EV*density/EPS0*np.clip(edges/scale,-a,a)
                    case['maximum_electric_error_over_boundary_field']=float(abs(e-exact).max()/abs(model.left_electric))
            history=dict(time_s=np.asarray(times),faces_m=np.asarray([s.faces for s in states]),
                         velocity_m_s=np.asarray([s.velocity for s in states]),thermal_j_m2_history=np.asarray([s.thermal for s in states]),**budgets)
            np.savez_compressed(out/(case['name']+'.npz'),numbers_m2=numbers,**history)
            report['cases'].append(case);finals[kind,cells]=(model,state,history);save()
            print('Separate supports',case['name'],'success',case['success'],'front error',case.get('maximum_front_error_m'),
                  'energy error',case['maximum_relative_energy_error'],flush=True)
    comparisons=[];grid=np.linspace(-.19,.19,20001)
    for low,high in zip(grids[:-1],grids[1:]):
        if not all(c['success'] for c in report['cases'] if c['kind']=='plasma' and c['cells'] in (low,high)):continue
        m,s,_=finals['plasma',low];n,t,_=finals['plasma',high]
        # J and E are piecewise linear on this union; two-point Gauss is exact
        # for their squared difference and includes the empty exterior.
        edges=np.unique(np.r_[m.bounds,n.bounds,s.faces.ravel(),t.faces.ravel()])
        nodes=np.array([.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)])
        points=(edges[:-1,None]+np.diff(edges)[:,None]*nodes).ravel()
        weights=np.repeat(np.diff(edges)/2,2)
        first=m.sample(s,points);second=n.sample(t,points)
        def relative(field):
            a,b=first[field],second[field];w=weights if a.ndim==1 else weights[:,None]
            return float(np.sqrt(np.sum(w*(a-b)**2)/np.sum(w*b*b)))
        comparisons.append(dict(coarse_cells=low,fine_cells=high,
            species_number_l1_over_initial=np.sum(weights*abs(first['number']-second['number']),axis=1).tolist(),
            electric_relative_l2=relative('electric'),current_relative_l2=relative('current'),
            species_front_difference_m=abs(s.faces[:,[0,-1]]-t.faces[:,[0,-1]]).max(axis=1).tolist()))
        comparisons[-1]['species_number_l1_over_initial']=(np.asarray(comparisons[-1]['species_number_l1_over_initial'])/(2*.01*1e14)).tolist()
    report['plasma_spatial_sensitivity']=comparisons
    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for cells in grids:
        for kind,ax in (('neutral',axes[0,0]),('charged',axes[0,1])):
            model,state,h=finals[kind,cells];x=.5*(state.faces[0,:-1]+state.faces[0,1:])
            density=1e19 if kind=='neutral' else 1e12
            ax.plot(x*1000,model.numbers[0]/np.diff(state.faces[0])/density,label=f'{cells} cells')
            if cells==grids[-1]:
                time_s=h['time_s'][-1];sound=np.sqrt((5/3)*EV*(.1 if kind=='neutral' else 1.)/model.mass[0])
                ax.plot(grid*1000,reference_density(grid,kind,time_s,.025,sound,1.5),'k--',label='Independent continuum reference')
        model,state,h=finals['plasma',cells];edges,e,_,_=model.field(state)
        completed=next(c['success'] for c in report['cases'] if c['name']==f'plasma_n{cells}')
        label=f'{cells} cells: {h["time_s"][-1]*1e9:.3f} ns'+('' if completed else ' STOPPED')
        axes[1,0].plot(edges*1000,e,label=label)
        axes[1,1].plot(h['time_s']*1e9,(h['energy_j_m2']-h['energy_j_m2'][0])/h['energy_j_m2'][0],label=label)
    for ax,title in zip(axes.flat,('Neutral vacuum rarefaction','Warm charged affine expansion','Electron-ion field','Total energy including empty exterior')):
        ax.set(title=title);ax.grid(alpha=.2);ax.legend(fontsize=7)
    for ax in axes[0]:ax.set(xlim=(-55,55),xlabel='Axial position [mm]',ylabel='Number / initial density')
    axes[1,0].set(xlim=(-90,90),xlabel='Axial position [mm]',ylabel='Electric field [V/m]')
    axes[1,1].set(xlabel='Time [ns]',ylabel='Energy change / initial energy')
    fig.suptitle('Independent species support: actual free-boundary fluid dynamics')
    fig.savefig(out/'species_support.png',dpi=170);fig.savefig(out/'species_support.pdf');plt.close(fig)
    if provenance(['examples.species_material'])!=identity:raise RuntimeError('Source changed during material species run')
    report['terminal']=True;report['success']=all(c['success'] for c in report['cases']);save();return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('--mesh',choices=('uniform','cosine'),default='uniform')
    p.add_argument('--viscosity',type=float,nargs=2,default=(0.,0.));p.add_argument('--step-scale',type=float,default=1.)
    p.add_argument('--viscosity-limiter',choices=('none','minmod'),default='none')
    p.add_argument('--cells',type=int,nargs='+',default=[32,64,128]);args=p.parse_args()
    result=run(args.destination,mesh=args.mesh,viscosity=args.viscosity,step_scale=args.step_scale,grids=args.cells,viscosity_limiter=args.viscosity_limiter);raise SystemExit(0 if result['success'] else 1)
