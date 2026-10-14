"""Actual moving-coil pulse histories, plasma fields and exterior reciprocity."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.material_cylinder import MaterialCylinder
from plasma_bench.material_two_species import MaterialTwoSpecies,MaterialSpeciesState
from plasma_bench.material_magnetic_vacuum import MaterialMagneticVacuum
from plasma_bench.material_vacuum_space import MaterialVacuumSpace
from plasma_bench.material_coil_current import MaterialCoilCurrent
from plasma_bench.material_carrier_current import MaterialCurrentState
from plasma_bench.material_vacuum_samples import physical_samples
from plasma_bench.material_plasma_samples import material_samples
from plasma_bench.waveforms import from_dict
from .implicit_plasma_records import sha


def reconstruct(path,d):
    mesh=MHDMesh(d['radius_m'],d['half_length_m'],d['nr'],d['nz']);q=d['quadrature'];p=d['degree']
    domain=MaterialCylinder(mesh,quadrature=q);space=MaterialVacuumSpace(mesh,degree=p,quadrature=q)
    r,z=space.reference_r,space.reference_z;profile=((1-(r/mesh.radius)**2)*(1-(z/mesh.half_length)**2)).reshape(d['nr']*d['nz'],-1)
    fluid=MaterialTwoSpecies(domain,domain.fluid_reference,d['central_density_m3']*profile**1.5,np.repeat(profile[...,None],2,axis=-1),d['ion_mass_kg'])
    fluid.initial.entropy[:]=-(2/3)*np.log(d['central_density_m3'])
    vacuum=MaterialMagneticVacuum(mesh,degree=p,quadrature_order=q,vacuum_cells=d['vacuum_cells'],modes=d['modes'])
    model=MaterialCoilCurrent(fluid,space,vacuum,d['coils'],from_dict(d['waveform']),resistivity=d['resistivity_ohm_m'])
    return fluid,model


def state_at(path,step):
    with np.load(path/('initial.npz' if step==0 else f'step_{step:05}.npz')) as raw:
        return MaterialCurrentState(MaterialSpeciesState(**{k:raw[k].copy() for k in MaterialSpeciesState.__dataclass_fields__}),raw['magnetic'].copy())


def build(destination,sources):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh physical plot directory required')
    runs=[]
    for source in sources:
        p=Path(source);d=json.loads((p/'summary.json').read_text())
        if not d['terminal'] or not d['success']:raise ValueError('Complete successful actual pulse required')
        runs.append((p,d))
    out.mkdir();fig,axes=plt.subplots(2,3,figsize=(15,9),constrained_layout=True);summary={}
    for p,d in runs:
        rows=[d['initial']]+d['records'];t=np.array([r['time'] for r in rows])*1e6;label=f"{d['nr']}×{d['nz']}, {d['steps']} steps, {d.get('integrator','rk4')}"
        axes[0,0].plot(t,[r['source_current_a'] for r in rows],label=label)
        j=np.array([r['peak_current_components_a_m2'] for r in rows]);axes[0,1].plot(t,j[:,1],label='Jθ '+label)
        axes[0,1].plot(t,np.max(j[:,[0,2]],axis=1),'--',label='Max Jr,Jz '+label)
        axes[1,0].plot(t,np.array([r['joule_heat_j'] for r in rows])*1e9,label='Joule heat '+label)
        axes[1,0].plot(t,np.array([r['magnetic_induced_energy_j'] for r in rows])*1e9,'--',label='Induced field '+label)
        axes[1,1].plot(t,[r['relative_external_work_residual'] for r in rows],label='External E work '+label)
        if all(r['relative_response_work_residual'] is not None for r in rows):
            axes[1,1].plot(t,[r['relative_response_work_residual'] for r in rows],'--',label='Mutual/source work '+label)
        rows=[d['initial']]+d['records'];summary[p.name]={'source_sha256':d['source_sha256'],'summary_sha256':sha(p/'summary.json'),
            'max_external_work_residual':max(abs(r['relative_external_work_residual']) for r in rows),
            'max_response_work_residual':max((abs(r['relative_response_work_residual']) for r in rows if r['relative_response_work_residual'] is not None),default=None),
            'integrator':d.get('integrator','rk4'),
            'max_current_components_a_m2':np.max([r['peak_current_components_a_m2'] for r in rows],axis=0).tolist(),
            'max_pair_power_residual_w':max(abs(r['pair_power_residual_w']) for r in rows),
            'minimum_entropy':np.min([r['minimum_specific_entropy'] for r in rows],axis=0).tolist(),'final':{k:v for k,v in rows[-1].items() if k!='boundary_coordinates'}}
    p,d=runs[-1];fluid,model=reconstruct(p,d);step=round((d['waveform']['rise']+.5*d['waveform']['hold'])/d['dt']);state=state_at(p,step)
    response=model.response(state,time=step*d['dt']);vol=response['fields']['volume'];coordinates=state.fluid.coordinates[fluid.domain.ids]
    j=E*response['fields']['n'].ravel()[:,None]*response['drift'];data=[response['magnetic_field'][:,2]*1e6,j[:,1]]
    for ax,value,title,unit in zip(axes[:,2],data,('Total axial field during the pulse','Induced azimuthal current during the pulse'),('Bz [µT]','Jθ [A/m²]')):
        color=(np.sum(vol*value.reshape(vol.shape),axis=1)/np.sum(vol,axis=1)).reshape(d['nr'],d['nz']);cap=float(abs(color).max())
        artist=ax.pcolormesh(coordinates[...,1]*100,coordinates[...,0]*100,color,cmap='RdBu_r',vmin=-cap,vmax=cap,rasterized=True)
        for zz,rr in ((coordinates[-1,:,1],coordinates[-1,:,0]),(coordinates[:,0,1],coordinates[:,0,0]),(coordinates[:,-1,1],coordinates[:,-1,0])):ax.plot(zz*100,rr*100,color='#233e4b',lw=1)
        for radius,center,turns in d['coils']:ax.plot(center*100,radius*100,'o',color='#aa6539',ms=7)
        ax.set(xlabel='Axial z [cm]',ylabel='Radius r [cm]',title=f"{title}\n{step*d['dt']*1e6:.3f} µs; copper dots: prescribed coils",aspect='equal');fig.colorbar(artist,ax=ax,shrink=.7,label=unit)
    coils=np.asarray(d['coils']);vac=response['vacuum_response'];psi,_,_=physical_samples(model.vacuum,state.fluid.coordinates,model.space.vacuum_trace@state.magnetic,coils[:,0],coils[:,1],geometry=vac)
    linkage=float(2*np.pi*np.sum(coils[:,2]*psi));volume_linkage=response['source_linkage_wb']
    summary['peak_pulse_exterior_reciprocity']=dict(time=step*d['dt'],source=p.name,volume_current_linkage_wb=volume_linkage,
        exterior_coil_flux_linkage_wb=linkage,relative_difference=linkage/volume_linkage-1,
        note='Independent physical reconstruction at coil positions; discrepancy is finite-space error, not a correction to source work.')
    if len(runs)==2 and all(runs[0][1][k]==runs[1][1][k] for k in ('nr','nz','degree','quadrature','duration','coils','waveform')):
        pairs=[reconstruct(path,record) for path,record in runs];comparisons=[]
        for fraction in (.25,.5,1.):
            samples=[];states=[]
            for (path,record),(f,m) in zip(runs,pairs):
                k=round(fraction*record['steps']);s=state_at(path,k);states.append(s);samples.append(m.response(s,rate=False,time=k*record['dt']))
            volume=samples[-1]['volume'];row={'time':fraction*d['duration']}
            for name,values in [('B',[v['magnetic_field'] for v in samples]),('J',[E*v['fields']['n'].ravel()[:,None]*v['drift'] for v in samples]),('velocity',[v['velocity'] for v in samples])]:
                row[name+'_relative_l2_difference']=float(np.sqrt(np.sum(volume[:,None]*(values[0]-values[1])**2)/np.sum(volume[:,None]*values[1]**2)))
                if name in ('B','J'):row[name+'_component_relative_l2_differences']=np.sqrt(np.sum(volume[:,None]*(values[0]-values[1])**2,axis=0)/np.sum(volume[:,None]*values[1]**2,axis=0)).tolist()
            row['maximum_position_difference_m']=float(abs(states[0].fluid.coordinates-states[1].fluid.coordinates).max());comparisons.append(row)
        same_method=runs[0][1].get('integrator','rk4')==runs[1][1].get('integrator','rk4')
        summary['same_grid_time_comparison' if same_method else 'same_grid_method_comparison']=comparisons
    elif len(runs)==2 and all(runs[0][1][k]==runs[1][1][k] for k in ('radius_m','half_length_m','central_density_m3','resistivity_ohm_m','duration','coils','waveform','ion_mass_kg')):
        pairs=[reconstruct(path,record) for path,record in runs];comparisons=[]
        radial=np.unique(np.concatenate([m.space.mesh.rf for f,m in pairs]))
        axial=np.unique(np.concatenate([m.space.mesh.zf for f,m in pairs]))
        x,w=np.polynomial.legendre.leggauss(10)
        r=(.5*(radial[1:]+radial[:-1])[:,None]+.5*np.diff(radial)[:,None]*x).ravel()
        z=(.5*(axial[1:]+axial[:-1])[:,None]+.5*np.diff(axial)[:,None]*x).ravel()
        rr,zz=np.meshgrid(r,z,indexing='ij')
        dv=2*np.pi*rr*(.5*np.diff(radial)[:,None]*w).ravel()[:,None]*(.5*np.diff(axial)[:,None]*w).ravel()[None,:]
        n0=d['central_density_m3']*((1-(rr/d['radius_m'])**2)*(1-(zz/d['half_length_m'])**2))**1.5
        for fraction in (.25,.5,1.):
            samples=[]
            for (path,record),(f,m) in zip(runs,pairs):
                k=round(fraction*record['steps']);s=state_at(path,k)
                samples.append(material_samples(m,s,rr,zz,n0,time=k*record['dt']))
            volume=dv*samples[-1]['volume_ratio'];row={'time':fraction*d['duration']}
            for name in ('magnetic','current','velocity'):
                a,b=[v[name] for v in samples]
                row[name+'_relative_l2_difference']=float(np.sqrt(np.sum(volume[...,None]*(a-b)**2)/np.sum(volume[...,None]*b*b)))
                if name!='velocity':row[name+'_component_relative_l2_differences']=np.sqrt(np.sum(volume[...,None]*(a-b)**2,axis=(0,1))/np.sum(volume[...,None]*b*b,axis=(0,1))).tolist()
            row['maximum_position_difference_m']=float(np.linalg.norm(samples[0]['points']-samples[1]['points'],axis=-1).max())
            comparisons.append(row)
        summary['material_label_spatial_comparison']=dict(reference_quadrature_points=rr.size,
            weighting='Ten-point Gaussian quadrature per union reference cell; actual volume of final listed run.',
            scope='Same material labels, each at its actual physical position; two grids do not establish convergence. Time errors remain.',states=comparisons)
    titles=['Prescribed finite-coil current','Actual current components','Heat and induced magnetic energy','Source-work balance']
    labels=['Common coil current [A]','Peak current density [A/m²]','Energy [nJ]','Work residual / initial plasma energy']
    for ax,title,label in zip(axes[:,:2].T.flatten(),[titles[0],titles[2],titles[1],titles[3]],[labels[0],labels[2],labels[1],labels[3]]):
        ax.set(xlabel='Time [µs]',ylabel=label,title=title);ax.grid(alpha=.2);ax.legend(fontsize=7,frameon=False)
    fig.suptitle('Finite coils drive a freely moving two-species plasma\nActual vacuum on side and ends · prescribed resistance · coarse, no formation claim',fontsize=15)
    fig.savefig(out/'material_coil_pulse.png',dpi=170);fig.savefig(out/'material_coil_pulse.pdf');plt.close(fig)
    (out/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n');print('moving coil figure complete',out,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('sources',nargs='+');a=p.parse_args();build(a.destination,a.sources)
