"""Read-only scientific fields and convergence from completed magnetic motion."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import mu_0 as MU0
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.material_two_species import MaterialTwoSpecies,MaterialSpeciesState
from plasma_bench.material_magnetic import CompactMagneticField,MaterialMagnetic,MaterialIdealPlasma,local_geometry
from .implicit_plasma_records import sha


def load(nr):
    path=Path(f'outputs/material_magnetic_r{nr}_t800_weak');record=json.loads((path/'summary.json').read_text())
    if not record['terminal'] or not record['success']:raise ValueError('Only terminal complete trajectories may be displayed')
    d=MovingDomain(nr,2*nr,radius=record['radius_m'],matching_radius=4*record['radius_m'],vacuum_layers=2,modes=1,quadrature=record['quadrature'])
    with np.load(path/'initial.npz') as saved:
        initial=MaterialSpeciesState(**{k:saved[k].copy() for k in MaterialSpeciesState.__dataclass_fields__});s=saved['profile_s']
    h=1-s*s;f=MaterialTwoSpecies(d,initial.coordinates,record['central_density_m3']*h**1.5,h[...,None]*record['central_temperature_ev'],record['ion_mass_kg'])
    m=MaterialMagnetic(f,CompactMagneticField(f,poloidal_tesla=record['poloidal_amplitude_t'],toroidal_tesla=record['toroidal_amplitude_t']))
    p=MaterialIdealPlasma(f,m)
    final_path=path/f"step_{record['steps']:05}.npz"
    with np.load(final_path) as saved:final=MaterialSpeciesState(**{k:saved[k].copy() for k in MaterialSpeciesState.__dataclass_fields__})
    return path,record,f,m,p,initial,final


def common(record,fluid,magnetic,state):
    nr=fluid.domain.nr;nt=fluid.domain.nt;sub=32//nr
    xi,tau=np.meshgrid((np.arange(sub)+.5)/sub,(np.arange(sub)+.5)/sub,indexing='ij');xi=xi.ravel();tau=tau.ravel()
    shape,x,gradient,det,inv=local_geometry(state.coordinates,fluid.ids,xi,tau)
    _,x0,_,det0,_=local_geometry(fluid.initial.coordinates,fluid.ids,xi,tau)
    b=magnetic._field_at(state,xi,tau)[0];f=fluid.fields(state)
    v=np.zeros_like(b);v[...,[0,2]]=np.einsum('qn,cna->cqa',shape,f['velocity'][fluid.ids]);v[...,1]=x[...,0]*np.einsum('qn,cn->cq',shape,f['omega'][fluid.ids])
    physical=magnetic.samples(state);j=np.einsum('qn,cna->cqa',shape,physical['current_coefficients'][fluid.ids]);j[...,:2]*=x[...,0,None]
    s=(np.repeat(np.arange(nr),nt)[:,None]+xi[None,:])/nr
    n0=record['central_density_m3']*(1-s*s)**1.5;ratio=det*x[...,0]/(det0*x0[...,0]);n=n0/ratio
    angle=np.einsum('qn,cn->cq',shape,state.angle[fluid.ids]);particles=n0*2*np.pi*x0[...,0]*det0/sub**2
    def order(a):return a.reshape(nr,nt,sub,sub,*a.shape[2:]).transpose(0,2,1,3,*range(4,4+len(a.shape[2:]))).reshape(32,64,*a.shape[2:])
    return {k:order(a) for k,a in dict(position=x,velocity=v,field=b,current=j,density=n,angle=angle,particles=particles).items()}


def build(destination):
    out=Path(destination)
    if out.exists():raise ValueError('Preserve existing figure outputs')
    data=[load(nr) for nr in (4,8,16)];out.mkdir(parents=True)
    comparison=[];common_data=[common(r,f,m,u) for _,r,f,m,p,i,u in data]
    for left,right in ((0,1),(1,2)):
        lo,hi=common_data[left],common_data[right];w=hi['particles'];w=w/w.sum();r=data[right][1]
        delta={}
        for key,scale in (('position',r['radius_m']),('velocity',r['speed_scale_m_s']),('field',r['poloidal_amplitude_t']),('current',r['poloidal_amplitude_t']/(MU0*r['radius_m'])),('density',r['central_density_m3']),('angle',1.)):
            error=lo[key]-hi[key];square=np.sum(error*error,axis=-1) if error.ndim==3 else error*error
            delta[key]=float(np.sqrt(np.sum(w*square))/scale)
        delta['grids']=[data[left][1]['radial_cells'],r['radial_cells']];comparison.append(delta)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(16,10),constrained_layout=True)
    path,r,f,m,p,initial,final=data[-1]
    current_limit=max(float(abs(m.samples(u)['current_density'][...,1]).max()) for u in (initial,final))/1e6
    for col,(state,label) in enumerate(((initial,'Initial'),(final,f"After {r['duration']*1e6:.2f} µs"))):
        g=p.fields(state);physical=m.samples(state);coords=state.coordinates[f.domain.ids[:f.domain.nr+1]]
        reference=initial.coordinates[f.domain.ids[:f.domain.nr+1]];s=np.arange(f.domain.nr+1)[:,None]/f.domain.nr
        psi=.5*r['poloidal_amplitude_t']*reference[...,0]**2*(1-s*s)**4
        for row,key,title,limit,cmap in ((0,'magnetic_field','Axial field Bz [T]',r['poloidal_amplitude_t'],'RdBu_r'),(1,'current_density','Toroidal current Jθ [MA/m²]',current_limit,'PuOr_r')):
            component=2 if row==0 else 1;z=physical[key][...,component]
            z=np.sum(z*g['volume'],axis=1)/np.sum(g['volume'],axis=1)
            if row:z=z/1e6
            z=z.reshape(f.domain.nr,f.domain.nt);limit=limit or max(abs(z).max(),1e-12)
            artist=axes[row,col].pcolormesh(coords[...,1]*100,coords[...,0]*100,z,cmap=cmap,vmin=-limit,vmax=limit,rasterized=True)
            axes[row,col].contour(coords[...,1]*100,coords[...,0]*100,psi,levels=np.linspace(psi.max()/8,psi.max()*7/8,7),colors='#132d36',linewidths=.7)
            axes[row,col].plot(coords[-1,:,1]*100,coords[-1,:,0]*100,color='#1a3842',linewidth=1.2)
            axes[row,col].set(xlabel='Axial z [cm]',ylabel='Radius r [cm]',title=f'{label} · {title}',aspect='equal',xlim=(-5.1,5.1),ylim=(0,5.2))
            fig.colorbar(artist,ax=axes[row,col],shrink=.72)
    colors=['#81939c','#b75f35','#167c80']
    for (path,r,f,m,p,i,u),color in zip(data,colors):
        history=[r['initial']]+r['records'];t=np.array([a['time'] for a in history])*1e6
        axes[0,2].plot(t,[a['magnetic_energy_j'] for a in history],label=f"{r['radial_cells']} × {r['angular_cells']}",color=color)
    axes[0,2].set(xlabel='Time [µs]',ylabel='Magnetic energy [J]',title='Magnetic energy converted into fluid motion');axes[0,2].legend(frameon=False);axes[0,2].grid(alpha=.2)
    meshes=np.array([d[1]['radial_cells'] for d in data])
    for key,label in (('initial_b_rms_error_over_amplitude','B / B₀'),('initial_current_rms_error_over_scale','J / (B₀/µ₀R)'),('initial_acceleration_particle_rms_error_over_scale','Acceleration / (C/R)')):
        axes[1,2].loglog(meshes,[d[1][key] for d in data],'o-',label=label)
    axes[1,2].set(xlabel='Radial cells',ylabel='Normalized RMS error',title='Independent spherical initial-field comparison');axes[1,2].legend(frameon=False);axes[1,2].grid(alpha=.2,which='both');axes[1,2].set_xticks(meshes,[str(v) for v in meshes])
    fig.suptitle('Moving magnetized fluid in exact vacuum\nIdeal-MHD coupling fixture · initialized closed flux · no coils, Hall response or resistive formation',fontsize=16)
    fig.savefig(out/'material_magnetic_expansion.png',dpi=180);fig.savefig(out/'material_magnetic_expansion.pdf');plt.close(fig)
    result=dict(material_label_rms_refinement=comparison,inputs={str(d[0]/'summary.json'):sha(d[0]/'summary.json') for d in data},
        interpretation='Complete ideal-MHD trajectories. Refinement norms compare the same material labels, not Eulerian positions. Current is weak Ampere including internal faces; no current inertia/Hall/coil formation claim.')
    (out/'comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':build('outputs/material_magnetic_comparison')
