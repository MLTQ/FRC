"""Read-only moving finite-mass current histories and physical final fields."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge as E,mu_0 as MU0,electron_mass as ME
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.material_cylinder import MaterialCylinder
from plasma_bench.material_current_space import MaterialCurrentSpace
from plasma_bench.material_two_species import MaterialTwoSpecies,MaterialSpeciesState
from plasma_bench.material_carrier_current import MaterialCarrierCurrent,MaterialCurrentState
from .implicit_plasma_records import sha
from plasma_bench.coulomb_transport import regime


def load(path):
    path=Path(path);record=json.loads((path/'summary.json').read_text())
    if not record['terminal'] or not record['success']:raise ValueError('Completed successful source run required')
    return path,record


def reconstruct(path,record):
    nr=record['nr'];nz=record['nz'];mesh=MHDMesh(record['radius_m'],record['half_length_m'],nr,nz)
    domain=MaterialCylinder(mesh,quadrature=6);space=MaterialCurrentSpace(mesh,degree=3,quadrature=6)
    r=space.reference_r.reshape(nr*nz,-1);z=space.reference_z.reshape(nr*nz,-1)
    h=(1-(r/mesh.radius)**2)*(1-(z/mesh.half_length)**2);n=record['central_density_m3']*h**1.5
    fluid=MaterialTwoSpecies(domain,domain.fluid_reference,n,np.broadcast_to(h[...,None],h.shape+(2,)),record['ion_mass_kg']);fluid.initial.entropy[:]=-(2/3)*np.log(record['central_density_m3'])
    model=MaterialCarrierCurrent(fluid,space);states=[]
    for name in ('initial.npz',f"step_{record['steps']:05}.npz"):
        with np.load(path/name) as raw:
            states.append(MaterialCurrentState(MaterialSpeciesState(**{k:raw[k].copy() for k in MaterialSpeciesState.__dataclass_fields__}),raw['magnetic'].copy()))
    return fluid,model,states



def compact_reference(r,z,d):
    R=d['radius_m'];L=d['half_length_m'];A=d['amplitude_t']
    hr=1-(r/R)**2;hz=1-(z/L)**2
    b=A*r/R*hr**3*hz**3
    j=np.stack((6*A*r*z/(MU0*R*L*L)*hr**3*hz**2,
        A/(MU0*R)*hz**3*(2*hr**3-6*(r/R)**2*hr**2)),axis=-1)
    n=d['central_density_m3']*(hr*hz)**1.5
    return b,j,n


def initial_comparison(path,d):
    fluid,model,states=reconstruct(path,d);a=model.response(states[0],rate=False)
    r=fluid.reference_geometry['points'][...,0].ravel();z=fluid.reference_geometry['points'][...,1].ravel()
    b,j,n=compact_reference(r,z,d);weight=a['volume'];current=E*n[:,None]*a['drift']
    result=dict(magnetic_relative_l2_error=float(np.sqrt(np.sum(weight*(a['magnetic_field']-b)**2)/np.sum(weight*b*b))),
        current_relative_l2_error=float(np.sqrt(np.sum(weight[:,None]*(current-j)**2)/np.sum(weight[:,None]*j*j))),
        drift_particle_relative_l2_error=float(np.sqrt(np.sum((weight*n)[:,None]*(a['drift']-j/(E*n[:,None]))**2)/np.sum((weight*n)[:,None]*(j/(E*n[:,None]))**2))))
    q,w=np.polynomial.legendre.leggauss(80);rv=.5*d['radius_m']*(q+1);zv=d['half_length_m']*q
    rr,zz=np.meshgrid(rv,zv,indexing='ij');b,j,n=compact_reference(rr,zz,d)
    dv=2*np.pi*rv[:,None]*(.5*d['radius_m']*w[:,None])*(d['half_length_m']*w[None,:])
    wb=float(np.sum(dv*b*b)/(2*MU0));mu=d['ion_mass_kg']*ME/(d['ion_mass_kg']+ME)
    kj=float(mu/(2*E*E)*np.sum(dv*np.sum(j*j,axis=-1)/n))
    result.update(reference_magnetic_energy_j=wb,reference_current_energy_j=kj,
        magnetic_energy_relative_error=a['magnetic_energy']/wb-1,current_energy_relative_error=a['current_energy']/kj-1)
    # Independent analytic current-stress divergence, retaining cylindrical
    # connection terms. Complex differentiation acts only on these formulas.
    r=fluid.reference_geometry['points'][...,0].ravel();z=fluid.reference_geometry['points'][...,1].ravel()
    b,j,n=compact_reference(r,z,d)
    def stress(rr,zz):
        bb,jj,nn=compact_reference(rr,zz,d)
        return mu/E**2*jj[..., :,None]*jj[...,None,:]/nn[...,None,None]
    epsilon=1e-25;tr=stress(r+1j*epsilon,z).imag/epsilon;tz=stress(r,z+1j*epsilon).imag/epsilon;value=stress(r,z)
    divergence=np.stack((tr[:,0,0]+tz[:,0,1]+value[:,0,0]/r,tr[:,1,0]+tz[:,1,1]+value[:,1,0]/r),axis=-1)
    hr=1-(r/d['radius_m'])**2;hz=1-(z/d['half_length_m'])**2;h=hr*hz
    pressure_force=10*E*d['central_density_m3']*h[:,None]**1.5*np.column_stack((r*hz/d['radius_m']**2,z*hr/d['half_length_m']**2))
    exact=(pressure_force+np.column_stack((-j[:,1]*b,j[:,0]*b))-divergence)/((d['ion_mass_kg']+ME)*n[:,None])
    rhs=model.derivative(states[0]);nodes=np.zeros_like(states[0].fluid.coordinates)
    nodes[fluid.active_r,0]=fluid.mass_factors[0].solve(rhs.fluid.momentum[fluid.active_r,0]);nodes[:,1]=fluid.mass_factors[1].solve(rhs.fluid.momentum[:,1])
    sampled=np.einsum('qn,cna->cqa',fluid.shape,nodes[fluid.ids]).reshape(-1,2)
    scale=10*E/((d['ion_mass_kg']+ME)*d['radius_m'])
    result['acceleration_particle_rms_error_over_C_R']=float(np.sqrt(np.sum(fluid.particles.ravel()*np.sum((sampled-exact)**2,axis=-1))/fluid.particles.sum())/scale)
    return result


def build(destination,free_paths):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh figure output required')
    refs=[load(f'outputs/material_current_reference_r{nr}_t400') for nr in (4,8)]
    free=[load(p) for p in free_paths];out.mkdir(parents=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(16,10),constrained_layout=True)
    for path,d in refs:
        rows=[d['initial']]+d['records'];t=np.array([r['time'] for r in rows])/d['duration']
        axes[0,0].plot(t,[r['peak_b_t']/d['initial']['peak_b_t'] for r in rows],label=f"Actual {d['nr']} × {d['nz']}")
        axes[1,0].semilogy(t,[r['magnetic_rms_error_over_initial_amplitude'] for r in rows],label=f"B: {d['nr']} × {d['nz']}")
        axes[1,0].semilogy(t,[r['current_rms_error_over_initial_scale'] for r in rows],'--',label=f"J: {d['nr']} × {d['nz']}")
    t=np.linspace(0,1,200);lam=1+t
    axes[0,0].plot(t,2/((1+lam)*lam**2),'k:',lw=2,label='Independent finite-inertia limit')
    axes[0,0].plot(t,lam**-2,color='#9e7771',ls='--',label='Ideal flux freezing')
    axes[0,0].set(xlabel='Time / expansion time',ylabel='Peak B / initial peak',title='Current inertia changes the evolving field');axes[0,0].legend(frameon=False,fontsize=8);axes[0,0].grid(alpha=.2)
    axes[1,0].set(xlabel='Time / expansion time',ylabel='RMS error / stated initial scale',title='Independent magnetic/current reference');axes[1,0].legend(frameon=False,fontsize=8);axes[1,0].grid(alpha=.2)
    for path,d in free:
        rows=[d['initial']]+d['records'];t=np.array([r['time'] for r in rows])*1e6
        axes[0,1].plot(t,np.array([r['magnetic_energy_j'] for r in rows])*1e9,label=f"Magnetic {d['nr']} × {d['nz']}")
        axes[0,1].plot(t,np.array([r['current_energy_j'] for r in rows])*1e9,'--',label=f"Current inertia {d['nr']} × {d['nz']}")
        axes[1,1].plot(t,[r['relative_total_energy_change'] for r in rows],label=f"{d['nr']} × {d['nz']}")
    axes[0,1].set(xlabel='Time [µs]',ylabel='Energy [nJ]',title='Actual free motion: magnetic and carrier energy');axes[0,1].legend(frameon=False,fontsize=8);axes[0,1].grid(alpha=.2)
    axes[1,1].set(xlabel='Time [µs]',ylabel='Relative total energy change',title='Unforced coupled-fluid energy');axes[1,1].legend(frameon=False);axes[1,1].grid(alpha=.2)
    path,d=free[-1];fluid,model,states=reconstruct(path,d);nr=d['nr'];nz=d['nz']
    responses=[model.response(u,rate=False) for u in states];limit=max(abs(v['magnetic_field']).max() for v in responses)*1e6
    regimes=[]
    for name,a in zip(('initial','final'),responses):
        density=a['fields']['n'].ravel();temps=a['fields']['temperature'].reshape(-1,2);drift=np.linalg.norm(a['drift'],axis=-1)
        for location,index in (('maximum_density',int(np.argmax(density))),('maximum_drift',int(np.argmax(drift))),('minimum_density',int(np.argmin(density)))):
            value=regime(density[index],*temps[index],d['ion_mass_kg'],field=a['magnetic_field'][index],current_density=E*density[index]*drift[index],length=d['radius_m'],timescale=d['duration'])
            value.update(state=name,sample=location,quadrature_index=index,reference_length_note='Body radius, not a resolved local gradient length; no local validity certificate.')
            if value['classical_coefficients'] is not None:
                collisions=value['classical_coefficients']['electron_ion_momentum_rate_s']*d['duration']
                value.update(electron_collisions_over_duration=collisions,zero_resistivity_review=bool(collisions>=.1))
            regimes.append(value)
    for row in range(2):
        state=states[-1];a=responses[-1];label=f"After {d['duration']*1e6:.2f} µs"
        coords=state.fluid.coordinates[fluid.domain.ids];vol=a['fields']['volume'];b=a['magnetic_field'].reshape(nr*nz,-1)
        cell_b=(np.sum(vol*b,axis=1)/np.sum(vol,axis=1)).reshape(nr,nz)*1e6
        current=E*a['fields']['n'][...,None]*a['drift'].reshape(nr*nz,-1,2)
        cell_jz=(np.sum(vol*current[...,1],axis=1)/np.sum(vol,axis=1)).reshape(nr,nz)/1e3
        color=cell_b if row==0 else cell_jz;cap=limit if row==0 else float(abs(cell_jz).max())
        artist=axes[row,2].pcolormesh(coords[...,1]*100,coords[...,0]*100,color,cmap='RdBu_r',vmin=-cap,vmax=cap,rasterized=True)
        x=np.sum(vol[...,None]*a['fields']['points'],axis=1)/np.sum(vol,axis=1)[:,None]
        current=E*a['fields']['n'][...,None]*a['drift'].reshape(nr*nz,-1,2);j=np.sum(vol[...,None]*current,axis=1)/np.sum(vol,axis=1)[:,None]
        # Arrow length is normalized only for visual direction; colors are Btheta.
        norm=np.linalg.norm(j,axis=-1);keep=norm>norm.max()*1e-3
        axes[row,2].quiver(x[keep,1]*100,x[keep,0]*100,j[keep,1]/norm[keep],j[keep,0]/norm[keep],color='#243b48',width=.004,scale=27)
        axes[row,2].plot(coords[-1,:,1]*100,coords[-1,:,0]*100,color='#203c49')
        axes[row,2].plot(coords[:,0,1]*100,coords[:,0,0]*100,color='#203c49');axes[row,2].plot(coords[:,-1,1]*100,coords[:,-1,0]*100,color='#203c49')
        initial_coords=states[0].fluid.coordinates[fluid.domain.ids]
        for zz,rr in ((initial_coords[-1,:,1],initial_coords[-1,:,0]),(initial_coords[:,0,1],initial_coords[:,0,0]),(initial_coords[:,-1,1],initial_coords[:,-1,0])):
            axes[row,2].plot(zz*100,rr*100,'k--',lw=.8,alpha=.5)
        component='Bθ [µT]' if row==0 else 'Jz [kA/m²]'
        axes[row,2].set(xlabel='Axial z [cm]',ylabel='Radius r [cm]',title=f'{label} · {component}; arrows: current direction',aspect='equal')
        fig.colorbar(artist,ax=axes[row,2],shrink=.7)
    fig.suptitle('Moving fluid with finite-mass carrier currents\nIsentropic toroidal-field sector · Hall/current stress included · no coils or FRC formation',fontsize=16)
    fig.savefig(out/'material_current_expansion.png',dpi=180);fig.savefig(out/'material_current_expansion.pdf');plt.close(fig)
    historical=[load(p) for p in ('outputs/material_current_free_r4_t2000','outputs/material_current_free_r8_t2000_refined')]
    summary={'historical_magnetic_only_initial_comparison':{p.name:initial_comparison(p,d) for p,d in historical},'free_initial_comparison':{p.name:initial_comparison(p,d) for p,d in free},'regime_samples':regimes,'inputs':{str(p/'summary.json'):sha(p/'summary.json') for p,d in refs+free+historical},'reference_final':{p.name:d['records'][-1] for p,d in refs},'free_final':{p.name:d['records'][-1] for p,d in free},'interpretation':'Reference motion is prescribed; free mode evolves simultaneous material/current coordinates. Constant entropy/eta=0 sector, no external coils, no formation claim.'}
    (out/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n');print('material current figure complete',out,flush=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('free_paths',nargs='+');a=p.parse_args();build(a.destination,a.free_paths)
