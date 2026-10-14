"""Read actual full-vector material states into an offline 3-D field/current view."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.constants import elementary_charge as E
from scipy.interpolate import LinearNDInterpolator
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.material_cylinder import MaterialCylinder
from plasma_bench.material_vector_space import MaterialVectorSpace
from plasma_bench.material_two_species import MaterialTwoSpecies,MaterialSpeciesState
from plasma_bench.material_vector_current import MaterialVectorCurrent
from plasma_bench.material_carrier_current import MaterialCurrentState
from .implicit_plasma_records import sha


def reconstruct(path,d):
    mesh=MHDMesh(d['radius_m'],d['half_length_m'],d['nr'],d['nz']);q=d['quadrature']
    domain=MaterialCylinder(mesh,quadrature=q);space=MaterialVectorSpace(mesh,degree=d['degree'],quadrature=q)
    r=space.reference_r.reshape(d['nr']*d['nz'],-1);z=space.reference_z.reshape(r.shape)
    profile=(1-(r/mesh.radius)**2)*(1-(z/mesh.half_length)**2)
    fluid=MaterialTwoSpecies(domain,domain.fluid_reference,d['central_density_m3']*profile**1.5,np.broadcast_to(profile[...,None],profile.shape+(2,)),d['ion_mass_kg'])
    fluid.initial.entropy[:]=-(2/3)*np.log(d['central_density_m3'])
    return fluid,MaterialVectorCurrent(fluid,space,resistivity=d['resistivity_ohm_m'],coulomb_exchange=d['coulomb_exchange'])


def trace(points,magnetic,radius,length):
    interpolate=LinearNDInterpolator(points,magnetic,fill_value=np.nan)
    seeds=np.array([[zz*length,rr*radius*np.cos(phi),rr*radius*np.sin(phi)] for zz in (-.45,0,.45) for rr in (.15,.35,.55,.75) for phi in np.arange(4)*np.pi/2])
    locations=seeds.copy();active=np.ones(len(seeds),dtype=bool);paths=[[] for _ in seeds];strength=[[] for _ in seeds]
    threshold=np.linalg.norm(magnetic,axis=-1).max()*1e-5;ds=radius/64
    def slope(x):
        r=np.hypot(x[:,1],x[:,2]);b=interpolate(np.column_stack((r,x[:,0])));norm=np.linalg.norm(b,axis=-1)
        valid=np.isfinite(b).all(axis=-1)&(norm>threshold)&(r>radius*1e-5);v=np.zeros_like(x)
        a=np.flatnonzero(valid)
        v[a]=np.column_stack((b[a,2],(b[a,0]*x[a,1]-b[a,1]*x[a,2])/r[a],(b[a,0]*x[a,2]+b[a,1]*x[a,1])/r[a]))/norm[a,None]
        return v,valid,norm
    for step in range(513):
        k1,m1,b=slope(locations);active&=m1
        for i in np.flatnonzero(active):paths[i].append(locations[i].tolist());strength[i].append(float(b[i]))
        if step==512 or not active.any():break
        k2,m2,_=slope(locations+.5*ds*k1);k3,m3,_=slope(locations+.5*ds*k2);k4,m4,_=slope(locations+ds*k3)
        active&=m2&m3&m4;locations[active]+=ds/6*(k1+2*k2+2*k3+k4)[active]
        if step>80:active&=np.linalg.norm(locations-seeds,axis=-1)>1.2*ds
    return [dict(points=p,field=b) for p,b in zip(paths,strength) if len(p)>5]


def outline(coordinates):
    paths=[];phi=np.linspace(0,2*np.pi,65)
    for r,z in coordinates[-1]:paths.append(np.column_stack((np.full_like(phi,z),r*np.cos(phi),r*np.sin(phi))))
    side=np.concatenate((coordinates[:,0],coordinates[-1,1:],coordinates[-2::-1,-1]))
    for angle in np.arange(8)*np.pi/4:paths.append(np.column_stack((side[:,1],side[:,0]*np.cos(angle),side[:,0]*np.sin(angle))))
    return paths


def build(source,destination):
    source=Path(source);out=Path(destination);d=json.loads((source/'summary.json').read_text())
    if out.exists() or not d['terminal'] or not d['success']:raise ValueError('Fresh display and a completed successful actual vector run required')
    fluid,model=reconstruct(source,d);frames=[];inputs={str(source/'summary.json'):sha(source/'summary.json')}
    steps=sorted(set((0,max(1,d['steps']//20),d['steps']//4,d['steps']//2,d['steps'])))
    q=d['quadrature'];local=np.array([q+1,2*q-2,(q-2)*q+1,(q-1)*q-2])
    indices=(np.arange(d['nr']*d['nz'])[:,None]*q*q+local[None,:]).ravel()
    for step in steps:
        filename='initial.npz' if step==0 else f'step_{step:05}.npz';inputs[str(source/filename)]=sha(source/filename)
        with np.load(source/filename) as raw:state=MaterialCurrentState(MaterialSpeciesState(**{k:raw[k].copy() for k in MaterialSpeciesState.__dataclass_fields__}),raw['magnetic'].copy())
        a=model.response(state,rate=False);points=a['fields']['points'].reshape(-1,2);current=E*a['fields']['n'].ravel()[:,None]*a['drift'];magnetic=a['magnetic_field']
        colors=np.column_stack((np.linalg.norm(current,axis=-1),current));selected=np.unique(np.r_[indices,np.argmax(abs(colors),axis=0)])
        cloud=[]
        for phi in np.arange(16)*np.pi/8:
            cloud.extend(np.column_stack((points[selected,1],points[selected,0]*np.cos(phi),points[selected,0]*np.sin(phi),colors[selected])).tolist())
        frames.append(dict(step=step,time_us=step*d['dt']*1e6,lines=trace(points,magnetic,d['radius_m'],d['half_length_m']),current=cloud,
            outline=[p.tolist() for p in outline(state.fluid.coordinates[fluid.domain.ids])],bmax=float(np.linalg.norm(magnetic,axis=-1).max()),jmax=abs(colors).max(axis=0).tolist(),
            temperature=(np.sum(fluid.particles[...,None]*a['fields']['temperature'],axis=(0,1))/fluid.particles.sum()).tolist()))
        print('vector display state',step,flush=True)
    extent=1.05*max(max(np.max(abs(np.asarray(p)[:,0])),np.max(abs(np.asarray(p)[:,1:]))/.6) for frame in frames for p in frame['outline']);physical_frames=frames
    out.mkdir(parents=True)
    fig=plt.figure(figsize=(14,7),constrained_layout=True);last=physical_frames[-1];axes=[fig.add_subplot(1,2,k+1,projection='3d') for k in range(2)]
    norm=Normalize(0,last['bmax']*1e6)
    for line in last['lines']:
        points=np.asarray(line['points'])*100;segments=np.stack((points[:-1],points[1:]),axis=1)
        artist=Line3DCollection(segments,cmap='viridis',norm=norm,linewidths=.8);artist.set_array(np.asarray(line['field'][1:])*1e6);axes[0].add_collection3d(artist)
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),ax=axes[0],shrink=.65,label='|B| [µT]')
    cloud=np.asarray(last['current']);cap=last['jmax'][0]
    artist=axes[1].scatter(cloud[:,0]*100,cloud[:,1]*100,cloud[:,2]*100,c=cloud[:,3],s=4,alpha=.65,cmap='viridis',vmin=0,vmax=cap)
    fig.colorbar(artist,ax=axes[1],shrink=.65,label='|J| [A/m²]')
    for ax,title in zip(axes,('Magnetic field lines','Current density magnitude')):
        for p in last['outline']:
            x=np.asarray(p)*100;ax.plot(*x.T,color='#6c7b82',lw=.35,alpha=.35)
        ax.set(xlabel='Axial z [cm]',ylabel='x [cm]',zlabel='y [cm]',title=title);ax.set_box_aspect((2,1.2,1.2));ax.view_init(22,-66)
        ax.set_xlim(-extent*100,extent*100);ax.set_ylim(-.6*extent*100,.6*extent*100);ax.set_zlim(-.6*extent*100,.6*extent*100)
    fig.suptitle(f"Actual full-vector state at {last['time_us']:.3f} µs\nInitialized closed flux · coarse {d['nr']}×{d['nz']} mesh · no coil-driven formation",fontsize=15)
    fig.savefig(out/'material_vector_3d.png',dpi=180);fig.savefig(out/'material_vector_3d.pdf');plt.close(fig)
    for frame in frames:
        for line in frame['lines']:line['points']=(np.asarray(line['points'])/extent).tolist()
        frame['outline']=[(np.asarray(p)/extent).tolist() for p in frame['outline']]
        cloud=np.asarray(frame['current']);cloud[:,:3]/=extent;frame['current']=cloud.tolist()
    data=dict(frames=frames,bmax=max(f['bmax'] for f in frames),jmax=np.max([f['jmax'] for f in frames],axis=0).tolist(),physical_scale_m=extent,
        scope=f"Source: {source.name}. Grid {d['nr']}×{d['nz']}, field degree {d['degree']}; {d['steps']} completed steps. Initial magnetic/current errors {100*d['initial_magnetic_relative_l2_error']:.1f}%/{100*d['initial_current_relative_l2_error']:.1f}%. Prescribed resistance, conditional Coulomb exchange and explicit upwind entropy diffusion. All B/J/flow components evolve; exterior field and matter are zero. Physical edge/transport validity and spatial convergence remain unresolved.")
    (out/'display_data.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    assets=Path(__file__).resolve().parents[1]/'visualizations/material_fields'
    html=(assets/'index.html').read_text().replace('/* DATA */','const DATA='+json.dumps(data,separators=(',',':'))+';').replace('/* VIEWER */',(assets/'viewer.js').read_text())
    (out/'index.html').write_text(html)
    (out/'source.json').write_text(json.dumps(dict(source=str(source),inputs=inputs,display_scope='Saved states only; linear spatial interpolation for paths, sampled current cloud, no temporal interpolation or formation inference.'),indent=2)+'\n')
    print('material vector display complete',out,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('destination');a=p.parse_args();build(a.source,a.destination)
