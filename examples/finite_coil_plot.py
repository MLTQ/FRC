"""Plot saved finite-coil trajectories and their actual current/field profiles."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.constants import electron_mass as ME
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.finite_cylinder_fluid import FiniteCylinderFluid
from plasma_bench.finite_cylinder_geometry import FiniteCylinderGeometry
from plasma_bench.finite_cylinder_map import FiniteCylinderMap
from plasma_bench.finite_coil_source import FiniteCoilSource
from plasma_bench.regular_current import RegularCurrent
from .current_transport_sampling import bases


def run(destination,coarse,fine):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=False)
    paths=[Path(coarse),Path(fine)];data=[json.loads((p/'summary.json').read_text()) for p in paths]
    if not all(d['terminal'] and d['success'] for d in data):raise ValueError('Only complete physical trajectories may be presented')
    d=data[-1];m=MHDMesh(d['radius_m'],d['length_m']/2,d['field_nr'],d['field_nz'])
    g=FiniteCylinderGeometry(m,quadrature_order=12,vacuum_cells=12,modes=8)
    fluid=FiniteCylinderFluid(m.radius,2*m.half_length,d['fluid_radial_elements'],d['fluid_axial_points'],d['ion_mass_kg']+ME)
    mapping=FiniteCylinderMap(fluid,g);source=FiniteCoilSource(g,d['coils'])
    row=min(d['records'],key=lambda x:abs(x['time']-6e-6));saved=np.load(paths[-1]/row['arrays_file']);p,toroidal=g.unpack(saved['magnetic'])
    current=RegularCurrent(mapping,saved['state'][...,0],d['ion_mass_kg']).current(saved['magnetic']);_,jt,_=g.split_electric(current)
    r=np.linspace(.0001,m.radius,120);z=np.linspace(-m.half_length,m.half_length,241);R=r[:,None];Z=z[None,:]
    nx,dx,_,_,_=bases(.5*m.rf**2,.5*r*r,5);nz,dz,_,_,_=bases(m.zf,z,5)
    psi_source,bsource=source.evaluate(R,Z)
    psi=nx@p@nz.T+float(saved['source_current'])*psi_source
    br=-(nx@p@dz.T)/R+float(saved['source_current'])*bsource[...,0]
    bz=dx@p@nz.T+float(saved['source_current'])*bsource[...,2]
    jtheta=(nx[:,1:]@(m.rf[1:,None]*jt.reshape(m.nr,m.nz+1))@nz.T)/R
    fig,axes=plt.subplots(2,2,figsize=(12,8.2),layout='constrained')
    for path,record,color in zip(paths,data,('#607d8b','#007f86')):
        rows=[record['initial'],*record['records']];times=np.array([x['time'] for x in rows])*1e6
        label=f"{record['field_nr']}×{record['field_nz']} field cells"
        axes[0,0].plot(times,np.array([x['mean_plasma_bz_t'] for x in rows])*1e6,color=color,label=label)
        axes[0,1].plot(times,np.array([x['maximum_current_a_m2'] for x in rows])/1e3,color=color,label=label)
    axes[0,0].plot(times,np.array([x['source_current_a'] for x in rows])*float(np.sum(g.weight*source.b_unit[...,2])/np.sum(g.weight))*1e6,'--',color='#cb6b28',label='Coil-only mean over same volume')
    axes[0,0].set(xlabel='Time (µs)',ylabel='Axial field (µT)',title='Field penetration — volume mean in plasma');axes[0,0].legend(fontsize=8)
    axes[0,1].set(xlabel='Time (µs)',ylabel='Peak |J| (kA/m²)',title='Induced plasma current');axes[0,1].legend(fontsize=8)
    vmax=float(abs(jtheta).max());pc=axes[1,0].pcolormesh(z,r*100,jtheta/1e3,cmap='RdBu_r',vmin=-vmax/1e3,vmax=vmax/1e3,shading='auto')
    fig.colorbar(pc,ax=axes[1,0],label='Jθ (kA/m²)');axes[1,0].set(xlabel='Axial position z (m)',ylabel='Radius (cm)',title=f'Actual current at {row["time"]*1e6:.1f} µs')
    pc=axes[1,1].pcolormesh(z,r*100,bz*1e6,cmap='viridis',shading='auto');fig.colorbar(pc,ax=axes[1,1],label='Bz (µT)')
    levels=np.linspace(float(psi.min()),float(psi.max()),20)[1:-1]
    axes[1,1].contour(z,r*100,psi,levels=levels,colors='white',linewidths=.6,alpha=.8)
    axes[1,1].set(xlabel='Axial position z (m)',ylabel='Radius (cm)',title='Actual Bz and poloidal magnetic paths')
    for ax in axes.ravel():ax.grid(alpha=.12)
    peak=d.get('peak_center_source_field_t',1e-4)*1e6
    fig.suptitle(f'Finite-coil plasma pulse • {peak:g} µT imposed center peak • evolving resistance',fontsize=15)
    fig.supxlabel('Sealed 1 m × 15 cm plasma cylinder; empty magnetic exterior. Axisymmetric response, not FRC formation.',fontsize=10)
    fig.savefig(destination/'finite_coil_response.png',dpi=180);fig.savefig(destination/'finite_coil_response.pdf')
    np.savez_compressed(destination/'profile.npz',r=r,z=z,current_theta=jtheta,br=br,bz=bz,psi=psi,time=row['time'])
    (destination/'inputs.json').write_text(json.dumps(dict(paths=[str(p) for p in paths],frame=row['arrays_file'],source_sha256=[d['source_sha256'] for d in data]),indent=2)+'\n')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('coarse');p.add_argument('fine');a=p.parse_args();run(a.destination,a.coarse,a.fine)
